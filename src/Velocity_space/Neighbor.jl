# Velocity-space neighbor index.
#
# A velocity grid (`AbstractVsData`) is a regular grid of `vs_trees_num[d]` root cells per
# dimension, each refined as an octree up to `AMR_VS_MAXLEVEL`, with every leaf carrying an
# explicit `level` and `midpoint`.  The leaves are stored in lexicographic-root /
# in-place-Morton order, which is *not* globally sorted, so a stand-alone index is built here
# to answer point-location and face-neighbor queries.
#
# Each leaf is keyed by `cellkey = root_lin * 2^(DIM*maxlevel) + intra_root_morton`, where
# `root_lin` is the row-major index of its root cell (a regular Cartesian lookup, so
# `vs_trees_num` may be *any* positive integers) and `intra_root_morton` is the Z-order code
# of its lower corner inside that root (the root subtree is a natural `2^maxlevel` octree).
# Within one root the octree leaves tile a `2^maxlevel`-per-dimension block, so each leaf owns
# a contiguous `cellkey` interval `[cellkey, cellkey + 2^(DIM*(maxlevel-level)))` that never
# spills across roots; the leaf containing a query point is therefore the predecessor of the
# point's `cellkey` (verified against that interval).  Build is `O(N log N)`, a query is
# `O(log N)`, storage is 8 B per velocity DOF, and the buffers are reusable across grids so the
# per-AMR-pass allocation is `O(max N_v)`.
#
# `vs_face_neighbor` is the reusable primitive: given a grid, a leaf and a face (`dim`, `dir`),
# it returns the same-or-coarser neighbor leaf across that face (`0` when the face is on the
# velocity-domain boundary, i.e. the vacuum tail).  It backs both the Löhner refinement
# indicator and, later, finite-difference velocity-space derivatives.

const _VS_INDEX_BITS = 24                      # low bits reserved for the leaf id
const _VS_INDEX_MASK = (UInt64(1) << _VS_INDEX_BITS) - 1

"""
$(TYPEDEF)
Reusable point-location / face-neighbor index over a single velocity grid.
$(TYPEDFIELDS)
"""
mutable struct VsNeighborIndex{DIM}
    """
    Lower corner of the velocity domain per dimension.
    """
    vmin::NTuple{DIM,Float64}
    """
    Finest-level cell size per dimension.
    """
    h_fine::NTuple{DIM,Float64}
    """
    Number of root cells per dimension (any positive integers).
    """
    vs_trees_num::NTuple{DIM,Int}
    """
    Row-major strides over the root grid (`rootstride[1] = 1`).
    """
    rootstride::NTuple{DIM,Int}
    """
    Maximum velocity-space refinement level.
    """
    maxlevel::Int
    """
    Bits of the intra-root Morton code, `DIM * maxlevel`.
    """
    mortonbits::Int
    """
    Sorted packed keys `(cellkey << INDEX_BITS) | leaf`; only the first `n` are valid.
    """
    keys::Vector{UInt64}
    """
    Reusable scratch buffer for the in-place sort (avoids per-build radix allocation).
    """
    scratch::Vector{UInt64}
    """
    Number of valid entries in `keys`.
    """
    n::Int
end

function VsNeighborIndex{DIM}() where {DIM}
    z = ntuple(_ -> 0, DIM)
    return VsNeighborIndex{DIM}(
        ntuple(_ -> 0.0, DIM),
        ntuple(_ -> 0.0, DIM),
        z,
        z,
        0,
        0,
        UInt64[],
        UInt64[],
        0,
    )
end

# Finest-level global integer coordinate of a leaf's lower corner from its center and level.
@inline function _corner_index(
    mid::Float64,
    level,
    vmin_d::Float64,
    h_fine_d::Float64,
    maxlevel::Int,
)
    cell = h_fine_d * (1 << (maxlevel - level))        # cell size along this dim
    return round(Int, (mid - 0.5 * cell - vmin_d) / h_fine_d)
end

"""
$(TYPEDSIGNATURES)
(Re)build `idx` over velocity grid `vs`.  `vmin`/`ds0` are the domain lower corner and the
root (level-0) cell size per dimension; `vs_trees_num` the per-dimension root count (any
positive integers) and `maxlevel` is `AMR_VS_MAXLEVEL`.  Reuses `idx`'s buffers when capacity
allows; the only work is one `O(N)` key fill plus the sort.
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
    vstn = ntuple(d -> Int(vs_trees_num[d]), Val(DIM))
    rootstride = ntuple(Val(DIM)) do d
        s = 1
        for e = 1:(d-1)
            s *= vstn[e]
        end
        s
    end
    mortonbits = DIM * maxlevel

    # Bit budget: cellkey (root_lin << mortonbits | morton) plus the leaf id must fit in 64 bits.
    prod_roots = 1
    @inbounds for d = 1:DIM
        prod_roots *= vstn[d]
    end
    rootbits = prod_roots <= 1 ? 0 : (sizeof(Int) * 8 - leading_zeros(prod_roots - 1))
    rootbits + mortonbits + _VS_INDEX_BITS > 64 && error(
        "VsNeighborIndex: cellkey ($(rootbits + mortonbits)) + index ($_VS_INDEX_BITS) bits exceed 64.",
    )
    N = vs.vs_num
    N > (1 << _VS_INDEX_BITS) - 1 &&
        error("VsNeighborIndex: $N leaves exceed the $_VS_INDEX_BITS-bit index field.")

    mask = (1 << maxlevel) - 1
    resize!(idx.keys, N)
    length(idx.scratch) < N && resize!(idx.scratch, N)
    keys = idx.keys
    @inbounds for i = 1:N
        L = vs.level[i]
        root_lin = 0
        morton = UInt64(0)
        for d = 1:DIM                     # inline; no tuple, no closure
            g = _corner_index(vs.midpoint[i, d], L, vmin[d], h_fine[d], maxlevel)
            root_lin += (g >> maxlevel) * rootstride[d]
            q = g & mask
            for b = 0:(maxlevel-1)
                morton |= (UInt64((q >> b) & 1)) << (b * DIM + (d - 1))
            end
        end
        cellkey = (UInt64(root_lin) << mortonbits) | morton
        keys[i] = (cellkey << _VS_INDEX_BITS) | UInt64(i)
    end
    sort!(
        keys;
        alg = Base.Sort.DEFAULT_UNSTABLE,
        order = Base.Order.Forward,
        scratch = idx.scratch,
    )

    idx.vmin = vmin
    idx.h_fine = h_fine
    idx.vs_trees_num = vstn
    idx.rootstride = rootstride
    idx.maxlevel = maxlevel
    idx.mortonbits = mortonbits
    idx.n = N
    return idx
end

# Locate the leaf owning `cellkey` (predecessor + interval containment check).  No allocation.
@inline function _vs_locate_code(
    idx::VsNeighborIndex{DIM},
    vs::AbstractVsData{DIM},
    code::UInt64,
) where {DIM}
    packed = (code << _VS_INDEX_BITS) | _VS_INDEX_MASK
    k = searchsortedlast(idx.keys, packed, 1, idx.n, Base.Order.Forward)
    k == 0 && return 0
    @inbounds begin
        corner = idx.keys[k] >> _VS_INDEX_BITS
        leaf = Int(idx.keys[k] & _VS_INDEX_MASK)
        L = vs.level[leaf]
        span = UInt64(1) << (DIM * (idx.maxlevel - L))   # owns [corner, corner+span), within its root
        (code >= corner && code < corner + span) && return leaf
    end
    return 0
end

"""
$(TYPEDSIGNATURES)
Return the leaf containing velocity-space `point`, or `0` if `point` lies outside the
velocity domain.  `vs` is the grid the index was built over (used to verify containment).
"""
function vs_locate(
    idx::VsNeighborIndex{DIM},
    vs::AbstractVsData{DIM},
    point::NTuple{DIM,Float64},
) where {DIM}
    mask = (1 << idx.maxlevel) - 1
    root_lin = 0
    morton = UInt64(0)
    @inbounds for d = 1:DIM
        g = floor(Int, (point[d] - idx.vmin[d]) / idx.h_fine[d])
        (g < 0 || g >= idx.vs_trees_num[d] << idx.maxlevel) && return 0
        root_lin += (g >> idx.maxlevel) * idx.rootstride[d]
        q = g & mask
        for b = 0:(idx.maxlevel-1)
            morton |= (UInt64((q >> b) & 1)) << (b * DIM + (d - 1))
        end
    end
    return _vs_locate_code(idx, vs, (UInt64(root_lin) << idx.mortonbits) | morton)
end

"""
$(TYPEDSIGNATURES)
Return the same-or-coarser neighbor leaf of leaf `i` across the face in dimension `dim`
and direction `dir` (`+1` / `-1`).  Returns `0` when the face is on the velocity-domain
boundary (treat as vacuum, `f = 0`).  Reusable primitive for Löhner indicators and
finite-difference velocity derivatives.
"""
function vs_face_neighbor(
    idx::VsNeighborIndex{DIM},
    vs::AbstractVsData{DIM},
    i::Integer,
    dim::Integer,
    dir::Integer,
) where {DIM}
    mask = (1 << idx.maxlevel) - 1
    root_lin = 0
    morton = UInt64(0)
    @inbounds begin
        L = vs.level[i]
        cell = idx.h_fine[dim] * (1 << (idx.maxlevel - L))
        # Probe just across the face, at the face-centre line in the other dimensions.
        for d = 1:DIM
            coord = vs.midpoint[i, d]
            d == dim && (coord += dir * (0.5 * cell + 0.5 * idx.h_fine[d]))
            g = floor(Int, (coord - idx.vmin[d]) / idx.h_fine[d])
            (g < 0 || g >= idx.vs_trees_num[d] << idx.maxlevel) && return 0
            root_lin += (g >> idx.maxlevel) * idx.rootstride[d]
            q = g & mask
            for b = 0:(idx.maxlevel-1)
                morton |= (UInt64((q >> b) & 1)) << (b * DIM + (d - 1))
            end
        end
    end
    return _vs_locate_code(idx, vs, (UInt64(root_lin) << idx.mortonbits) | morton)
end
