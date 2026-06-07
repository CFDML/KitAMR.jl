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
    "Number of valid entries in `keys`."
    n::Int
end

function VsNeighborIndex{DIM}() where {DIM}
    return VsNeighborIndex{DIM}(
        ntuple(_ -> 0.0, DIM), ntuple(_ -> 0.0, DIM), 0, 0, 0, UInt64[], 0,
    )
end

# Interleave DIM integer coordinates (nbits each) into a Morton code.
@inline function _morton_encode(q::NTuple{DIM,Int}, nbits::Int) where {DIM}
    code = UInt64(0)
    @inbounds for b in 0:nbits-1
        for d in 1:DIM
            code |= (UInt64((q[d] >> b) & 1)) << (b * DIM + (d - 1))
        end
    end
    return code
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
root (level-0) cell size per dimension; `maxlevel` is `AMR_VS_MAXLEVEL`.  Reuses `idx`'s
buffers when capacity allows.
"""
function build_vs_index!(
    idx::VsNeighborIndex{DIM},
    vs::AbstractVsData{DIM,NDF},
    vmin::NTuple{DIM,Float64},
    ds0::NTuple{DIM,Float64},
    maxlevel::Integer,
) where {DIM,NDF}
    maxlevel = Int(maxlevel)
    h_fine = ntuple(d -> ds0[d] / 2.0^maxlevel, DIM)
    N = vs.vs_num

    # Finest-grid size per dimension: octree leaves are aligned, so the largest covered
    # corner index gives the grid extent.  Require it to be a power of two and cubic.
    nfine_d = fill(0, DIM)
    @inbounds for i in 1:N
        L = vs.level[i]
        span = 1 << (maxlevel - L)
        for d in 1:DIM
            top = _corner_index(vs.midpoint[i, d], L, vmin[d], h_fine[d], maxlevel) + span
            top > nfine_d[d] && (nfine_d[d] = top)
        end
    end
    nbits = -1
    @inbounds for d in 1:DIM
        _ispow2(nfine_d[d]) ||
            error("VsNeighborIndex: finest grid size $(nfine_d[d]) along dim $d is not a power of two (need power-of-two vs_trees_num).")
        b = trailing_zeros(nfine_d[d])
        if nbits == -1
            nbits = b
        elseif nbits != b
            error("VsNeighborIndex: non-cubic velocity grid $(nfine_d); equal bits per dimension required.")
        end
    end
    nbits * DIM + _VS_INDEX_BITS > 64 &&
        error("VsNeighborIndex: morton ($(nbits*DIM)) + index ($_VS_INDEX_BITS) bits exceed 64.")
    N > (1 << _VS_INDEX_BITS) - 1 &&
        error("VsNeighborIndex: $N leaves exceed the $_VS_INDEX_BITS-bit index field.")

    length(idx.keys) < N && resize!(idx.keys, N)
    keys = idx.keys
    @inbounds for i in 1:N
        L = vs.level[i]
        q = ntuple(d -> _corner_index(vs.midpoint[i, d], L, vmin[d], h_fine[d], maxlevel), DIM)
        keys[i] = (_morton_encode(q, nbits) << _VS_INDEX_BITS) | UInt64(i)
    end
    sort!(view(keys, 1:N))

    idx.vmin = vmin
    idx.h_fine = h_fine
    idx.nbits = nbits
    idx.maxlevel = maxlevel
    idx.ncoarse = 1 << nbits
    idx.n = N
    return idx
end

"""
$(TYPEDSIGNATURES)
Return the leaf containing velocity-space `point`, or `0` if `point` lies outside the
velocity domain.  `vs` is the grid the index was built over (used to verify containment).
"""
function vs_locate(idx::VsNeighborIndex{DIM}, vs::AbstractVsData{DIM}, point::NTuple{DIM,Float64}) where {DIM}
    q = ntuple(d -> floor(Int, (point[d] - idx.vmin[d]) / idx.h_fine[d]), DIM)
    @inbounds for d in 1:DIM
        (q[d] < 0 || q[d] >= idx.ncoarse) && return 0
    end
    target = _morton_encode(q, idx.nbits)
    packed = (target << _VS_INDEX_BITS) | _VS_INDEX_MASK
    k = searchsortedlast(view(idx.keys, 1:idx.n), packed)
    k == 0 && return 0
    @inbounds begin
        corner = idx.keys[k] >> _VS_INDEX_BITS
        leaf = Int(idx.keys[k] & _VS_INDEX_MASK)
        L = vs.level[leaf]
        span = UInt64(1) << (DIM * (idx.maxlevel - L))
        # The predecessor owns [corner, corner+span); confirm the point falls inside.
        (target >= corner && target < corner + span) && return leaf
    end
    return 0
end

"""
$(TYPEDSIGNATURES)
Return the same-or-coarser neighbor leaf of leaf `i` across the face in dimension `dim`
and direction `dir` (`+1` / `-1`).  Returns `0` when the face is on the velocity-domain
boundary (treat as vacuum, `f = 0`).  Reusable primitive for Löhner indicators and
finite-difference velocity derivatives.
"""
function vs_face_neighbor(idx::VsNeighborIndex{DIM}, vs::AbstractVsData{DIM}, i::Integer, dim::Integer, dir::Integer) where {DIM}
    @inbounds begin
        L = vs.level[i]
        cell = idx.h_fine[dim] * (1 << (idx.maxlevel - L))
        # Probe just across the face, at the face-centre line in the other dimensions.
        point = ntuple(d -> d == dim ?
                            vs.midpoint[i, d] + dir * (0.5 * cell + 0.5 * idx.h_fine[d]) :
                            vs.midpoint[i, d], DIM)
    end
    return vs_locate(idx, vs, point)
end
