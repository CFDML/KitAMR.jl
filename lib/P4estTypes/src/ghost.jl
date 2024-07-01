"""
    GhostLayer{X,P}

Stores a ghost layer of quadrants that neighbor the domain local to the rank
for a `Pxest{X}`.  Also stores the corresponding local domain quadrants,
mirrors, that are in other rank's ghost layers.

# Fields
$(DocStringExtensions.FIELDS)

# See also
- [`ghostlayer`](@ref): a function used to construct a `GhostLayer`
"""
mutable struct GhostLayer{X,P}
    """The pointer (of type `P`) can be a pointer to either a
    `P4estTypes.P4est.p4est_ghost_t` or a
    `P4estTypes.P4est.p8est_ghost_t`.  See the help
    documentation for these types for more information about the
    underlying p4est structures. """
    pointer::P
    function GhostLayer{4}(pointer::Ptr{p4est_ghost_t})
        ghost = new{4,typeof(pointer)}(pointer)
        finalizer(ghost) do p
            p4est_ghost_destroy(p.pointer)
            p.pointer = C_NULL
            return
        end
    end
    function GhostLayer{8}(pointer::Ptr{p8est_ghost_t})
        ghost = new{8,typeof(pointer)}(pointer)
        finalizer(ghost) do p
            p8est_ghost_destroy(p.pointer)
            p.pointer = C_NULL
            return
        end
    end
end

"""
    ghostlayer(forest::Pxest{X}; connection=P4estTypes.CONNECT_FULL(Val(X)))

Construct a ghost layer of quadrants that neighbor the local to the rank
for the given `forest`.  Here `connection` determines what neighboring
quadrants to include (across face, edge, corner, or full) and can take the
values:
- `P4estTypes.CONNECT_FULL(Val(4))`: get face and corner neighbors.
- `P4estTypes.CONNECT_FULL(Val(8))`: get face, edge, and corner neighbors.
- `P4estTypes.CONNECT_FACE(Val(4))`: get face neighbors.
- `P4estTypes.CONNECT_FACE(Val(8))`: get face neighbors.
- `P4estTypes.CONNECT_EDGE(Val(8))`: get face and edge neighbors.
- `P4estTypes.CONNECT_CORNER(Val(4))`: get face and corner neighbors.
- `P4estTypes.CONNECT_CORNER(Val(8)): `get face, edge, and corner neighbors.
"""
function ghostlayer(forest::Pxest{X}; connection = CONNECT_FULL(Val(X))) where {X}
    return GhostLayer{X}((pxest_ghost_new(Val(X)))(forest, connection))
end

struct Ghosts{X,P,Q} <: AbstractArray{QuadrantWrapper,1}
    pointer::P
    ghostlayer::GhostLayer{X,Q}
end

Base.size(g::Ghosts) = (PointerWrapper(g.pointer).elem_count[],)
Base.@propagate_inbounds function Base.getindex(g::Ghosts{X}, i::Int) where {X}
    @boundscheck checkbounds(g, i)
    GC.@preserve g begin
        Q = pxest_quadrant_t(Val(X))
        quadrant = Ptr{Q}(pointer(PointerWrapper(g.pointer).array) + sizeof(Q) * (i - 1))
        return QuadrantWrapper{X,Ptr{Q}}(quadrant)
    end
end
Base.IndexStyle(::Ghosts) = IndexLinear()

struct Mirrors{X,P,Q} <: AbstractArray{QuadrantWrapper,1}
    pointer::P
    ghostlayer::GhostLayer{X,Q}
end

Base.size(m::Mirrors) = (PointerWrapper(m.pointer).elem_count[],)
Base.@propagate_inbounds function Base.getindex(m::Mirrors{X}, i::Int) where {X}
    @boundscheck checkbounds(m, i)
    GC.@preserve m begin
        Q = pxest_quadrant_t(Val(X))
        quadrant = Ptr{Q}(pointer(PointerWrapper(m.pointer).array) + sizeof(Q) * (i - 1))
        return QuadrantWrapper{X,Ptr{Q}}(quadrant)
    end
end
Base.IndexStyle(::Mirrors) = IndexLinear()

"""
    ghosts(gl::GhostLayer)

Returns an array-like structure with the [`QuadrantWrapper`](@ref)s that neighbor the
domain of the local rank.
"""
function ghosts(gl::GhostLayer{X}) where {X}
    gs = pointer(PointerWrapper(gl.pointer).ghosts)
    return Ghosts{X,typeof(gs),typeof(gl.pointer)}(gs, gl)
end

"""
    mirrors(gl::GhostLayer)

Returns an array-like structure with the [`QuadrantWrapper`](@ref)s in the local
domain that are in neighboring rank's ghost layers.
"""
function mirrors(gl::GhostLayer{X}) where {X}
    ms = pointer(PointerWrapper(gl.pointer).mirrors)
    return Mirrors{X,typeof(ms),typeof(gl.pointer)}(ms, gl)
end

function Base.unsafe_convert(
    ::Type{Ptr{p4est_ghost_t}},
    p::GhostLayer{4,Ptr{p4est_ghost_t}},
)
    return p.pointer
end

function Base.unsafe_convert(
    ::Type{Ptr{p8est_ghost_t}},
    p::GhostLayer{8,Ptr{p8est_ghost_t}},
)
    return p.pointer
end

"""
    expand!(ghost, forest, nodes)

Expand the ghost layer to include all elements with nodes shared with
neighboring ranks.

Consider the following forest-of-quadtrees

    +-------+---+---+
    |       | 5 | 6 |
    |   2   +---+---+
    |       | 3 | 4 |
    +-------+---+---+
    |       |       |
    |   0   |   1   |
    |       |       |
    +-------+-------+

that is partitioned so rank 0 owns quadrants {0,1} and rank 1 owns quadrants
{2, 3, 4, 5, 6}.  A fully connected ghost layer on rank 0 would include
quadrants {2, 3, 4}.  QuadrantWrapper 5 shares a global node with rank 0 but is not
in the fully connected ghost layer.  This function expands the ghost layer
to include quadrants like this.

See [`sharers`](@ref) to get a list of the global nodes shared with neighboring
ranks.
"""
function expand!(ghost::GhostLayer{X}, forest::Pxest{X}, nodes::LNodes{X}) where {X}
    pxest_ghost_support_lnodes(Val(X))(forest, nodes, ghost)
end

"""
    unsafe_mirror_proc_offsets(ghost::GhostLayer)

Returns 0-based indices into `mirror_proc_mirrors` for each rank.  This includes
an extra entry at the end of the array so that 1-based range into
`mirror_proc_mirrors` for rank `r` can be built with
```
(mirror_proc_offsets[r]+1):mirror_proc_offsets[r+1]
```

See `@doc P4estTypes.P4est.p4est_ghost_t` and
`@doc P4estTypes.P4est.p8est_ghost_t` for a more details.

Note, this unsafely wraps a C array.  So, you must ensure that the `ghost`
structure is preserved while using the return value.
"""
function unsafe_mirror_proc_offsets(ghost::GhostLayer{X}) where {X}
    gp = PointerWrapper(ghost.pointer)
    ptr = pointer(gp.mirror_proc_offsets)
    sz = (Int(gp.mpisize[] + 1),)
    mirror_proc_offsets = unsafe_wrap(Vector{eltype(ptr)}, ptr, sz, own = false)
    return mirror_proc_offsets
end

"""
    unsafe_proc_offsets(ghost::GhostLayer)

Returns 0-based indices into [`ghosts`](@ref) for each rank.  This includes
an extra entry at the end of the array so that 1-based range into
[`ghosts`](@ref) for rank `r` can be built with
```
(proc_offsets[r]+1):proc_offsets[r+1]
```
Thus the ghost quadrants associated with rank `r` can be obtained with
```
ghosts(ghost)[(proc_offsets[r]+1):proc_offsets[r+1]]
```

See `@doc P4estTypes.P4est.p4est_ghost_t` and
`@doc P4estTypes.P4est.p8est_ghost_t` for a more details.

Note, this unsafely wraps a C array.  So, you must ensure that the `ghost`
structure is preserved while using the return value.
"""
function unsafe_proc_offsets(ghost::GhostLayer{X}) where {X}
    gp = PointerWrapper(ghost.pointer)
    ptr = pointer(gp.proc_offsets)
    sz = (Int(gp.mpisize[] + 1),)
    proc_offsets = unsafe_wrap(Vector{eltype(ptr)}, ptr, sz, own = false)
    return proc_offsets
end

"""
    unsafe_mirror_proc_mirrors(ghost::GhostLayer)

Returns 0-based indices into [`mirrors`](@ref).  This is used in conjunction
with `mirror_proc_offsets` to get the mirror quadrants associated with
each rank.  For example
```
rrange = (mirror_proc_offsets[r]+1):mirror_proc_offsets[r+1]
mirrors(ghost)[mirror_proc_mirrors(rrange)]
```
selects the mirror quadrants associated with rank `r`.

See `@doc P4estTypes.P4est.p4est_ghost_t` and
`@doc P4estTypes.P4est.p8est_ghost_t` for a more details.

Note, this unsafely wraps a C array.  So, you must ensure that the `ghost`
structure is preserved while using the return value.
"""
function unsafe_mirror_proc_mirrors(ghost::GhostLayer{X}) where {X}
    gp = PointerWrapper(ghost.pointer)
    ptr = pointer(gp.mirror_proc_mirrors)
    sz = (Int(gp.mirror_proc_offsets[gp.mpisize[]+1]),)
    mirror_proc_mirrors = unsafe_wrap(Vector{eltype(ptr)}, ptr, sz, own = false)
    return mirror_proc_mirrors
end
