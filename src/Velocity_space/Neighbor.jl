# Velocity-space neighbor index.
#
# A velocity grid (`AbstractVsData`) is a linear forest of octrees: `vs_trees_num[d]`
# root cells per dimension, each refined as a `2^DIM`-tree, with every leaf carrying an
# explicit `level` and `midpoint`.  The leaves are stored in lexicographic-root /
# in-place-Morton order, which is *not* globally sorted, so a stand-alone index is built
# here to answer point-location and face-neighbor queries.
#
# The index is a single sorted `UInt64` array packing `(morton_code << INDEX_BITS) | leaf`,
# where `morton_code` is the Z-order code of the leaf's lower corner at the finest level.
# Because octree leaves tile the domain, the leaf owning a query point is the predecessor
# of the point's Morton code (verified geometrically).  Build is `O(N log N)`, a query is
# `O(log N)`, and storage is 8 B per velocity DOF.  The buffers are reusable across grids
# (`build_vs_index!`) so the per-AMR-pass allocation is `O(max N_v)`, independent of the
# total grid size.
#
# `vs_face_neighbor` is the reusable primitive: given a grid, a leaf and a face
# (`dim`, `dir`), it returns the same-or-coarser neighbor leaf across that face (`0` when
# the face is on the velocity-domain boundary, i.e. the vacuum tail).  It is intended to
# back both the Löhner refinement indicator and, later, finite-difference velocity-space
# derivatives.

const _VS_INDEX_BITS = 24                      # low bits reserved for the leaf id
const _VS_INDEX_MASK = (UInt64(1) << _VS_INDEX_BITS) - 1

"""
$(TYPEDEF)
Reusable point-location / face-neighbor index over a single velocity grid.
$(TYPEDFIELDS)
"""
mutable struct VsNeighborIndex{DIM}
    "Lower corner of the velocity domain per dimension."
    vmin::NTuple{DIM,Float64}
    "Finest-level cell size per dimension."
    h_fine::NTuple{DIM,Float64}
    "Bits per dimension, `log2(vs_trees_num[d] * 2^maxlevel)` (assumed equal across dims)."
    nbits::Int
    "Maximum velocity-space refinement level."
    maxlevel::Int
    "Number of finest cells per dimension, `2^nbits`."
    ncoarse::Int
    "Sorted packed keys `(morton << INDEX_BITS) | leaf`; only the first `n` are valid."
    keys::Vector{UInt64}
    "Reusable scratch buffer for the in-place sort (avoids per-build radix allocation)."
    scratch::Vector{UInt64}
    "Number of valid entries in `keys`."
    n::Int
end

function VsNeighborIndex{DIM}() where {DIM}
    return VsNeighborIndex{DIM}(
        ntuple(_ -> 0.0, DIM), ntuple(_ -> 0.0, DIM), 0, 0, 0, UInt64[], UInt64[], 0,
    )
end

@inline _ispow2(x::Integer) = x > 0 && (x & (x - 1)) == 0

# Finest-level integer coordinate of a leaf's lower corner from its center and level.
@inline function _corner_index(mid::Float64, level, vmin_d::Float64, h_fine_d::Float64, maxlevel::Int)
    cell = h_fine_d * (1 << (maxlevel - level))        # cell size along this dim
    return round(Int, (mid - 0.5 * cell - vmin_d) / h_fine_d)
end

"""
$(TYPEDSIGNATURES)
(Re)build `idx` over velocity grid `vs`.  `vmin`/`ds0` are the domain lower corner and the
root (level-0) cell size per dimension; `vs_trees_num` the per-dimension root count and
`maxlevel` is `AMR_VS_MAXLEVEL`.  The finest-grid bit width is derived directly from
`vs_trees_num`/`maxlevel` (no scan).  Reuses `idx`'s key buffer when capacity allows; the
only work is one `O(N)` key fill plus the sort.
"""
function build_vs_index!(
    idx::VsNeighborIndex{DIM},
    vs::AbstractVsData{DIM,NDF},
    vmin::NTuple{DIM,Float64},
    ds0::NTuple{DIM,Float64},
    vs_trees_num,
    maxlevel::Integer,
) where {DIM,NDF}
    maxlevel = Int(maxlevel)
    h_fine = ntuple(d -> ds0[d] / 2.0^maxlevel, Val(DIM))

    # Bits per dimension = log2(vs_trees_num[d]) + maxlevel; require power-of-two and cubic.
    nbits = -1
    @inbounds for d in 1:DIM
        nt = Int(vs_trees_num[d])
        _ispow2(nt) ||
            error("VsNeighborIndex: vs_trees_num[$d] = $nt is not a power of two.")
        b = trailing_zeros(nt) + maxlevel
        if nbits == -1
            nbits = b
        elseif nbits != b
            error("VsNeighborIndex: non-cubic velocity grid (vs_trees_num = $vs_trees_num); equal bits per dimension required.")
        end
    end
    nbits * DIM + _VS_INDEX_BITS > 64 &&
        error("VsNeighborIndex: morton ($(nbits*DIM)) + index ($_VS_INDEX_BITS) bits exceed 64.")
    N = vs.vs_num
    N > (1 << _VS_INDEX_BITS) - 1 &&
        error("VsNeighborIndex: $N leaves exceed the $_VS_INDEX_BITS-bit index field.")

    resize!(idx.keys, N)
    length(idx.scratch) < N && resize!(idx.scratch, N)
    keys = idx.keys
    @inbounds for i in 1:N
        L = vs.level[i]
        code = UInt64(0)
        for d in 1:DIM                     # inline Morton accumulation; no tuple, no closure
            qd = _corner_index(vs.midpoint[i, d], L, vmin[d], h_fine[d], maxlevel)
            for b in 0:nbits-1
                code |= (UInt64((qd >> b) & 1)) << (b * DIM + (d - 1))
            end
        end
        keys[i] = (code << _VS_INDEX_BITS) | UInt64(i)
    end
    sort!(keys; alg = Base.Sort.DEFAULT_UNSTABLE, order = Base.Order.Forward,
          scratch = idx.scratch)

    idx.vmin = vmin
    idx.h_fine = h_fine
    idx.nbits = nbits
    idx.maxlevel = maxlevel
    idx.ncoarse = 1 << nbits
    idx.n = N
    return idx
end

# Locate the leaf owning Morton code `code` (predecessor + geometric containment check).
# No view / tuple / closure allocation.
@inline function _vs_locate_code(idx::VsNeighborIndex{DIM}, vs::AbstractVsData{DIM}, code::UInt64) where {DIM}
    packed = (code << _VS_INDEX_BITS) | _VS_INDEX_MASK
    k = searchsortedlast(idx.keys, packed, 1, idx.n, Base.Order.Forward)
    k == 0 && return 0
    @inbounds begin
        corner = idx.keys[k] >> _VS_INDEX_BITS
        leaf = Int(idx.keys[k] & _VS_INDEX_MASK)
        L = vs.level[leaf]
        span = UInt64(1) << (DIM * (idx.maxlevel - L))   # predecessor owns [corner, corner+span)
        (code >= corner && code < corner + span) && return leaf
    end
    return 0
end

"""
$(TYPEDSIGNATURES)
Return the leaf containing velocity-space `point`, or `0` if `point` lies outside the
velocity domain.  `vs` is the grid the index was built over (used to verify containment).
"""
function vs_locate(idx::VsNeighborIndex{DIM}, vs::AbstractVsData{DIM}, point::NTuple{DIM,Float64}) where {DIM}
    code = UInt64(0)
    @inbounds for d in 1:DIM
        qd = floor(Int, (point[d] - idx.vmin[d]) / idx.h_fine[d])
        (qd < 0 || qd >= idx.ncoarse) && return 0
        for b in 0:idx.nbits-1
            code |= (UInt64((qd >> b) & 1)) << (b * DIM + (d - 1))
        end
    end
    return _vs_locate_code(idx, vs, code)
end

"""
$(TYPEDSIGNATURES)
Return the same-or-coarser neighbor leaf of leaf `i` across the face in dimension `dim`
and direction `dir` (`+1` / `-1`).  Returns `0` when the face is on the velocity-domain
boundary (treat as vacuum, `f = 0`).  Reusable primitive for Löhner indicators and
finite-difference velocity derivatives.
"""
function vs_face_neighbor(idx::VsNeighborIndex{DIM}, vs::AbstractVsData{DIM}, i::Integer, dim::Integer, dir::Integer) where {DIM}
    code = UInt64(0)
    @inbounds begin
        L = vs.level[i]
        cell = idx.h_fine[dim] * (1 << (idx.maxlevel - L))
        # Probe just across the face, at the face-centre line in the other dimensions.
        for d in 1:DIM
            coord = vs.midpoint[i, d]
            d == dim && (coord += dir * (0.5 * cell + 0.5 * idx.h_fine[d]))
            qd = floor(Int, (coord - idx.vmin[d]) / idx.h_fine[d])
            (qd < 0 || qd >= idx.ncoarse) && return 0
            for b in 0:idx.nbits-1
                code |= (UInt64((qd >> b) & 1)) << (b * DIM + (d - 1))
            end
        end
    end
    return _vs_locate_code(idx, vs, code)
end
