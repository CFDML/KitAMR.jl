@inline pxest_t(::Val{4}) = p4est_t
@inline pxest_t(::Val{8}) = p8est_t
@inline pxest_new_ext(::Val{4}) = p4est_new_ext
@inline pxest_new_ext(::Val{8}) = p8est_new_ext
@inline pxest_balance_ext(::Val{4}) = p4est_balance_ext
@inline pxest_balance_ext(::Val{8}) = p8est_balance_ext
@inline pxest_coarsen_ext(::Val{4}) = p4est_coarsen_ext
@inline pxest_coarsen_ext(::Val{8}) = p8est_coarsen_ext
@inline pxest_refine_ext(::Val{4}) = p4est_refine_ext
@inline pxest_refine_ext(::Val{8}) = p8est_refine_ext
@inline pxest_partition_ext(::Val{4}) = p4est_partition_ext
@inline pxest_partition_ext(::Val{8}) = p8est_partition_ext
@inline pxest_partition_lnodes(::Val{4}) = p4est_partition_lnodes
@inline pxest_partition_lnodes(::Val{8}) = p8est_partition_lnodes
@inline pxest_ghost_t(::Val{4}) = p4est_ghost_t
@inline pxest_ghost_t(::Val{8}) = p8est_ghost_t
@inline pxest_ghost_new(::Val{4}) = p4est_ghost_new
@inline pxest_ghost_new(::Val{8}) = p8est_ghost_new
@inline pxest_lnodes_t(::Val{4}) = p4est_lnodes
@inline pxest_lnodes_t(::Val{8}) = p8est_lnodes
@inline pxest_lnodes_new(::Val{4}) = p4est_lnodes_new
@inline pxest_lnodes_new(::Val{8}) = p8est_lnodes_new
@inline pxest_quadrant_t(::Val{4}) = p4est_quadrant
@inline pxest_quadrant_t(::Val{8}) = p8est_quadrant
@inline pxest_tree_t(::Val{4}) = p4est_tree
@inline pxest_tree_t(::Val{8}) = p8est_tree
@inline pxest_iter_volume_info_t(::Val{4}) = p4est_iter_volume_info_t
@inline pxest_iter_volume_info_t(::Val{8}) = p8est_iter_volume_info_t
@inline pxest_ghost_support_lnodes(::Val{4}) = p4est_ghost_support_lnodes
@inline pxest_ghost_support_lnodes(::Val{8}) = p8est_ghost_support_lnodes

# TODO move vtk wrapping functions to P4est.jl
function _p4est_vtk_context_new(forest, prefix)
    @ccall P4est.LibP4est.libp4est.p4est_vtk_context_new(
        forest::Ptr{p4est_t},
        prefix::Cstring,
    )::Ptr{Cvoid}
end
function _p8est_vtk_context_new(forest, prefix)
    @ccall P4est.LibP4est.libp4est.p8est_vtk_context_new(
        forest::Ptr{p8est_t},
        prefix::Cstring,
    )::Ptr{Cvoid}
end
@inline pxest_vtk_context_new(::Val{4}) = _p4est_vtk_context_new
@inline pxest_vtk_context_new(::Val{8}) = _p8est_vtk_context_new

function _p4est_vtk_context_set_geom(context, geometry)
    @ccall P4est.LibP4est.libp4est.p4est_vtk_context_set_geom(
        context::Ptr{Cvoid},
        geometry::Ptr{Cvoid},
    )::Cvoid
end
function _p8est_vtk_context_set_geom(context, geometry)
    @ccall P4est.LibP4est.libp4est.p8est_vtk_context_set_geom(
        context::Ptr{Cvoid},
        geometry::Ptr{Cvoid},
    )::Cvoid
end
@inline pxest_vtk_context_set_geom(::Val{4}) = _p4est_vtk_context_set_geom
@inline pxest_vtk_context_set_geom(::Val{8}) = _p8est_vtk_context_set_geom

function _p4est_vtk_context_set_scale(context, scale)
    @ccall P4est.LibP4est.libp4est.p4est_vtk_context_set_scale(
        context::Ptr{Cvoid},
        scale::Cdouble,
    )::Cvoid
end
function _p8est_vtk_context_set_scale(context, scale)
    @ccall P4est.LibP4est.libp4est.p8est_vtk_context_set_scale(
        context::Ptr{Cvoid},
        scale::Cdouble,
    )::Cvoid
end
@inline pxest_vtk_context_set_scale(::Val{4}) = _p4est_vtk_context_set_scale
@inline pxest_vtk_context_set_scale(::Val{8}) = _p8est_vtk_context_set_scale

function _p4est_vtk_context_set_continuous(context, continuous)
    @ccall P4est.LibP4est.libp4est.p4est_vtk_context_set_continuous(
        context::Ptr{Cvoid},
        continuous::Cint,
    )::Cvoid
end
function _p8est_vtk_context_set_continuous(context, continuous)
    @ccall P4est.LibP4est.libp4est.p8est_vtk_context_set_continuous(
        context::Ptr{Cvoid},
        continuous::Cint,
    )::Cvoid
end
@inline pxest_vtk_context_set_continuous(::Val{4}) = _p4est_vtk_context_set_continuous
@inline pxest_vtk_context_set_continuous(::Val{8}) = _p8est_vtk_context_set_continuous

function _p4est_vtk_write_header(context)
    @ccall P4est.LibP4est.libp4est.p4est_vtk_write_header(context::Ptr{Cvoid})::Ptr{Cvoid}
end
function _p8est_vtk_write_header(context)
    @ccall P4est.LibP4est.libp4est.p8est_vtk_write_header(context::Ptr{Cvoid})::Ptr{Cvoid}
end
@inline pxest_vtk_write_header(::Val{4}) = _p4est_vtk_write_header
@inline pxest_vtk_write_header(::Val{8}) = _p8est_vtk_write_header

function _p4est_vtk_write_footer(context)
    @ccall P4est.LibP4est.libp4est.p4est_vtk_write_footer(context::Ptr{Cvoid})::Cint
end
function _p8est_vtk_write_footer(context)
    @ccall P4est.LibP4est.libp4est.p8est_vtk_write_footer(context::Ptr{Cvoid})::Cint
end
@inline pxest_vtk_write_footer(::Val{4}) = _p4est_vtk_write_footer
@inline pxest_vtk_write_footer(::Val{8}) = _p8est_vtk_write_footer

function _p4est_vtk_write_cell_data(
    context,
    write_tree,
    write_level,
    write_rank,
    wrap_rank,
    num_cell_scalars,
    num_cell_vectors,
    fieldnames,
    values,
)
    @ccall P4est.LibP4est.libp4est.p4est_vtk_write_cell_data(
        context::Ptr{Cvoid},
        write_tree::Cint,
        write_level::Cint,
        write_rank::Cint,
        wrap_rank::Cint,
        num_cell_scalars::Cint,
        num_cell_vectors::Cint,
        fieldnames::Ptr{Cvoid},
        values::Ptr{Cvoid},
    )::Ptr{Cvoid}
end
function _p8est_vtk_write_cell_data(
    context,
    write_tree,
    write_level,
    write_rank,
    wrap_rank,
    num_cell_scalars,
    num_cell_vectors,
    fieldnames,
    values,
)
    @ccall P4est.LibP4est.libp4est.p8est_vtk_write_cell_data(
        context::Ptr{Cvoid},
        write_tree::Cint,
        write_level::Cint,
        write_rank::Cint,
        wrap_rank::Cint,
        num_cell_scalars::Cint,
        num_cell_vectors::Cint,
        fieldnames::Ptr{Cvoid},
        values::Ptr{Cvoid},
    )::Ptr{Cvoid}
end
@inline pxest_vtk_write_cell_data(::Val{4}) = _p4est_vtk_write_cell_data
@inline pxest_vtk_write_cell_data(::Val{8}) = _p8est_vtk_write_cell_data

"""
    P4estTypes.CONNECT_FULL(::Val{4})

Returns an integer indicating connecting quadrants across faces and corners.
"""
@inline CONNECT_FULL(::Val{4}) = P4EST_CONNECT_FULL
"""
    P4estTypes.CONNECT_FULL(::Val{8})

Returns an integer indicating connecting octants across faces, edges,
and corners.
"""
@inline CONNECT_FULL(::Val{8}) = P8EST_CONNECT_FULL
"""
    P4estTypes.CONNECT_FACE(::Val{4})

Returns an integer indicating connecting quadrants across faces.
"""
@inline CONNECT_FACE(::Val{4}) = P4EST_CONNECT_FACE
"""
    P4estTypes.CONNECT_FACE(::Val{8})

Returns an integer indicating connecting octants across faces.
"""
@inline CONNECT_FACE(::Val{8}) = P8EST_CONNECT_FACE
"""
    P4estTypes.CONNECT_CORNER(::Val{4})

Returns an integer indicating connecting quadrants across faces and corners.
"""
@inline CONNECT_CORNER(::Val{4}) = P4EST_CONNECT_CORNER
"""
    P4estTypes.CONNECT_CORNER(::Val{8})

Returns an integer indicating connecting octants across faces, edges,
and corners.
"""
@inline CONNECT_CORNER(::Val{8}) = P8EST_CONNECT_CORNER
"""
    P4estTypes.CONNECT_EDGE(::Val{8})

Returns an integer indicating connecting octants across faces and edges.
"""
@inline CONNECT_EDGE(::Val{8}) = P8EST_CONNECT_EDGE

const Locidx = P4est.p4est_locidx_t
const Gloidx = P4est.p4est_gloidx_t

"""
    QuadrantWrapper{X,P}

Stores a Pxest{X} quadrant (where `X=4` indicates a quadrant
and `X=8` indicates an octant; quadrant is used both as the general
term and the term for the 2D object).

# Fields

$(DocStringExtensions.FIELDS)
"""
struct QuadrantWrapper{X,P}
    """The pointer (of type `P`) can be a pointer to either a
   `P4estTypes.P4est.p4est_quadrant` or a
   `P4estTypes.P4est.p8est_quadrant`.  See the help
   documentation for these types for more information about the
   underlying p4est structures. """
    pointer::P
end

"""
    level(quadrant::QuadrantWrapper)

Returns the level of refinement for the quadrant.  Level 0 is the coarsest
level and `P4estTypes.P4est.P4EST_QMAXLEVEL` is the maximum refinement level.
"""
@inline level(quadrant::QuadrantWrapper) =
    GC.@preserve quadrant PointerWrapper(quadrant.pointer).level[]

"""
    coordinates(quadrant::QuadrantWrapper{4})

Returns a tuple of the quadrant's integer coordinates inside its tree.
"""
@inline function coordinates(quadrant::QuadrantWrapper{4})
    GC.@preserve quadrant begin
        qs = PointerWrapper(quadrant.pointer)
        return (qs.x[], qs.y[])
    end
end

"""
    unsafe_which_tree(quadrant::QuadrantWrapper)

Returns the `which_tree` field of the underlying quadrant.  This value is only
sometimes set so the function is marked unsafe.
"""
@inline function unsafe_which_tree(quadrant::QuadrantWrapper)
    return GC.@preserve quadrant PointerWrapper(quadrant.pointer).p.piggy3.which_tree[] +
                                 0x1
end

"""
    coordinates(quadrant::QuadrantWrapper{8})

Returns a tuple of the quadrant's integer coordinates inside its tree.
"""
@inline function coordinates(quadrant::QuadrantWrapper{8})
    GC.@preserve quadrant begin
        qs = PointerWrapper(quadrant.pointer)
        return (qs.x[], qs.y[], qs.z[])
    end
end

"""
    unsafe_storeuserdata!(quadrant::QuadrantWrapper, data)

Store the user data `data` associated with the `quadrant`.
"""
function unsafe_storeuserdata!(quadrant::QuadrantWrapper{X}, data::T) where {X,T}
    GC.@preserve quadrant begin
        qs = PointerWrapper(quadrant.pointer)
        unsafe_store!(Ptr{T}(qs.p.user_data[]), data)
    end
end

"""
    unsafe_loaduserdata(quadrant::QuadrantWrapper, type::Type)

Return the user data of type `type` associated with the `quadrant`.
"""
function unsafe_loaduserdata(quadrant::QuadrantWrapper{X}, ::Type{T}) where {X,T}
    GC.@preserve quadrant begin
        unsafe_load(Ptr{T}(PointerWrapper(quadrant.pointer).p.user_data[]))
    end
end

"""
    unsafe_local_num(quadrant::QuadrantWrapper)

Returns the `local_num` field of the underlying quadrant.  This value is only
sometimes set so the function is marked unsafe.  For example, it is set for
quadrants returned by [`ghosts`](@ref) and [`mirrors`](@ref).
"""
@inline function unsafe_local_num(quadrant::QuadrantWrapper)
    return GC.@preserve quadrant PointerWrapper(quadrant.pointer).p.piggy3.local_num[] + 0x1
end

"""
    Tree{X,P,Q} <: AbstractArray{QuadrantWrapper,1}

Stores the quadrants in a tree of a Pxest{X}.

# Fields

$(DocStringExtensions.FIELDS)
"""
struct Tree{X,P,Q} <: AbstractArray{QuadrantWrapper,1}
    """The pointer (of type `P`) can be a pointer to either a
    `P4estTypes.P4est.p4est_tree` or a
    `P4estTypes.P4est.p8est_tree`.  See the help documentation
    for these types for more information about the underlying
    p4est structures."""
    pointer::P
    """The forest (of type `Q`) the tree is associated with.  This is stored
    so the forest will not be reclaimed by the garbage collector too early.
    """
    forest::Q
end

Base.size(t::Tree) =
    (convert(Int, (GC.@preserve t PointerWrapper(t.pointer).quadrants.elem_count[])),)
function Base.getindex(t::Tree{X}, i::Int) where {X}
    @boundscheck checkbounds(t, i)
    GC.@preserve t begin
        Q = pxest_quadrant_t(Val(X))
        quadrant =
            Ptr{Q}(pointer(PointerWrapper(t.pointer).quadrants.array) + sizeof(Q) * (i - 1))
        return QuadrantWrapper{X,Ptr{Q}}(quadrant)
    end
end
Base.IndexStyle(::Tree) = IndexLinear()

"""
    offset(tree::Tree)

The cumulative sum of the quadrants over earlier trees on this rank (locals
only).
"""
offset(tree::Tree) = GC.@preserve tree PointerWrapper(tree.pointer).quadrants_offset[]

"""
    Pxest{X,P,C} <: AbstractArray{P4estTypes.Tree,1}

Stores the forest of quadtrees (when `X=4`) or octrees (when `X=8`).

This forest of octrees can be accessed in two ways.  First, as an array-of-arrays.
Each rank holds an array of quadrants for each tree of the [`Connectivity`](@ref)
associated with the forest. (Note, the quadrants are distributed among the ranks.
So, each rank will only have access to the quadrants it owns.) Second, using
`iterateforest` to iterate over the volumes, faces, edges, and corners of the
forest via callback functions.

# Fields

$(DocStringExtensions.FIELDS)

# See also

  - [`pxest`](@ref): a function that constructs a `Pxest` from a [`Connectivity`](@ref).
  - [`iterateforest`](@ref): a function to iterate over the volumes, faces, edges, and
    corners of the forest.
  - [`refine!`](@ref): refine the quadrants of the forest.
  - [`coarsen!`](@ref): coarsen the quadrants of the forest.
  - [`balance!`](@ref): two-to-one balance the quadrants of the forest.
  - [`partition!`](@ref): partition the quadrants of the forest.
  - [`ghostlayer`](@ref): get the ghost layer of quadrants for the forest.
  - [`lnodes`](@ref): get a global node numbering.
  - [`P4estTypes.savevtk`](@ref): save a VTK representation of the forest.
"""
mutable struct Pxest{X,P,C} <: AbstractArray{Tree,1}
    """The pointer (of type `P`) can be a pointer to either a
    `P4estTypes.P4est.LibP4est.p4est` or a
    `P4estTypes.P4est.LibP4est.p8est`.  See the help documentation for these
    types for more information about the underlying p4est structures."""
    pointer::P
    """The connectivity (of type `C`) the forest is associated with.  This is
    stored so the connectivity will not be reclaimed by the garbage collector
    too early."""
    connectivity::C
    """The MPI Communicator that includes the ranks participating in the
    forest."""
    comm::MPI.Comm
    function Pxest{4}(pointer::Ptr{p4est_t}, connectivity::Connectivity{4}, comm::MPI.Comm)
        forest = new{4,typeof(pointer),typeof(connectivity)}(pointer, connectivity, comm)
        finalizer(forest) do p
            p4est_destroy(p.pointer)
            p.pointer = C_NULL
            return
        end
    end
    function Pxest{8}(pointer::Ptr{p8est_t}, connectivity::Connectivity{8}, comm::MPI.Comm)
        forest = new{8,typeof(pointer),typeof(connectivity)}(pointer, connectivity, comm)
        finalizer(forest) do p
            p8est_destroy(p.pointer)
            p.pointer = C_NULL
            return
        end
    end
end

"""
    pxest(connectivity::Connectivity{X}; kw...) where {X}

Generate a distributed forest of quadtrees (if `X=4`) or octrees (if `X=8`)
based on `connectivity`.  Each element of `connectivity` becomes a tree root.

The connectivity is duplicated on all ranks but the leaves of the forest
are split (based on a space-filling curve order) among the ranks.

The keyword arguments (`kw...`) that control the construction of the forest
are:

  - `comm = MPI.COMM_WORLD`: the MPI Communicator object of the ranks sharing
    the forest.
  - `min_quadrants = 0`: the minimum number of quadrants per rank.  (This makes
    the initial refinement pattern `MPI.Comm_size` specific.)
  - `min_level = 0`: the minimum level of quadrant refinement for the forest.
  - `fill_uniform = true`: if `true` the forest will be filled with a uniform
    mesh otherwise it is the coarsest possible mesh.
  - `data_type = Nothing`: an `isbitstype` of the user data stored for each
    quadrant.
  - `init_function = nothing`: callback function with
    prototype `init_function(forest, treeid, quadrant)` called for each quadrant
    to initialized the user data.
"""
function pxest(
    connectivity::Connectivity{X};
    comm = MPI.COMM_WORLD,
    min_quadrants = 0,
    min_level = 0,
    fill_uniform = true,
    data_type = Nothing,
    init_function = nothing,
) where {X}
    MPI.Initialized() || MPI.Init()

    pointer = (pxest_new_ext(Val(X)))(
        comm,
        connectivity,
        min_quadrants,
        min_level,
        fill_uniform,
        sizeof(data_type),
        C_NULL,
        C_NULL,
    )

    forest = Pxest{X}(pointer, connectivity, comm)

    if !isnothing(init_function)
        init(forest, _, quadrant, _, treeid, _) = init_function(forest, treeid, quadrant)
        iterateforest(forest; volume = init)
    end

    return forest
end

"""
    quadrantstyle(forest)

Returns 4 if the `forest` quadrants are 2D and 8 if they are 3D.
"""
quadrantstyle(::Pxest{X}) where {X} = X

"""
    quadrantndims(forest)

Returns 2 if the `forest` quadrants are 2D and 3 if they are 3D.
"""
function quadrantstyle end
quadrantndims(::Pxest{4}) = 2
quadrantndims(::Pxest{8}) = 3

"""
    lengthoflocalquadrants(forest)

Return the number of quadrants local to the rank.
"""
lengthoflocalquadrants(p::Pxest) =
    GC.@preserve p PointerWrapper(p.pointer).local_num_quadrants[]
"""
    lengthofglobalquadrants(forest)

Return the number of quadrants distributed across the whole forest.
"""
lengthofglobalquadrants(p::Pxest) =
    GC.@preserve p PointerWrapper(p.pointer).global_num_quadrants[]
"""
    comm(forest)

Return the MPI Communicator used by the forest.
"""
comm(p::Pxest) = p.comm
"""
    connectivity(forest)

Return the [`Connectivity`](@ref) used by the forest.
"""
connectivity(p::Pxest) = p.connectivity

function Base.unsafe_convert(::Type{Ptr{p4est_t}}, p::Pxest{4,Ptr{p4est_t}})
    return p.pointer
end
function Base.unsafe_convert(::Type{Ptr{p8est_t}}, p::Pxest{8,Ptr{p8est_t}})
    return p.pointer
end

Base.size(p::Pxest) =
    (convert(Int, (GC.@preserve p PointerWrapper(p.pointer).trees.elem_count[])),)
function Base.getindex(p::Pxest{X}, i::Int) where {X}
    @boundscheck checkbounds(p, i)
    GC.@preserve p begin
        TR = pxest_tree_t(Val(X))
        tree =
            Ptr{TR}(pointer(PointerWrapper(p.pointer).trees.array) + sizeof(TR) * (i - 1))
        return Tree{X,Ptr{TR},typeof(p)}(tree, p)
    end
end
Base.IndexStyle(::Pxest) = IndexLinear()

"""
    unsafe_global_first_quadrant(forest::Pxest)

Returns 0-based indices into global quadrants.  This includes
an extra entry at the end of the array so that 1-based range into
the global quadrants for rank `r` can be built with

```
(global_first_quadrant[r]+1):global_first_quadrant[r+1]
```

Note, this unsafely wraps a C array.  So, you must ensure that the `forest`
structure is preserved while using the return value.
"""
@inline function unsafe_global_first_quadrant(forest::Pxest)
    fp = PointerWrapper(forest.pointer)

    return unsafe_wrap(
        Vector{p4est_gloidx_t},
        pointer(fp.global_first_quadrant),
        (fp.mpisize[] + 0x1,),
        own = false,
    )
end

function iterate_volume_callback(info, _)
    info = PointerWrapper(info)
    data = unsafe_pointer_to_objref(pointer(info.p4est.user_pointer))[]
    X = quadrantstyle(data.forest)
    quadrant = QuadrantWrapper{X,Ptr{pxest_quadrant_t(Val(X))}}(pointer(info.quad))
    data.volume(
        data.forest,
        data.ghost,
        quadrant,
        info.quadid[] + 0x1,
        info.treeid[] + 0x1,
        data.userdata,
    )
    return
end

@generated function generate_volume_callback(::Val{X}) where {X}
    I = pxest_iter_volume_info_t(Val(X))
    quote
        @cfunction(iterate_volume_callback, Cvoid, (Ptr{$I}, Ptr{Cvoid}))
    end
end

function generate_volume_closure(::Val{4}, volume, forest, ghost, userdata)
    function c(info::Ptr{p4est_iter_volume_info_t}, ::Ptr{Cvoid})::Cvoid
        info = PointerWrapper(info)

        quadrant = QuadrantWrapper{4,Ptr{p4est_quadrant_t}}(pointer(info.quad))

        volume(forest, ghost, quadrant, info.quadid[] + 0x1, info.treeid[] + 0x1, userdata)

        return
    end
    @cfunction($c, Cvoid, (Ptr{p4est_iter_volume_info_t}, Ptr{Cvoid}))
end

function generate_volume_closure(::Val{8}, volume, forest, ghost, userdata)
    function c(info::Ptr{p8est_iter_volume_info_t}, ::Ptr{Cvoid})::Cvoid
        info = PointerWrapper(info)

        quadrant = QuadrantWrapper{8,Ptr{p8est_quadrant_t}}(pointer(info.quad))

        volume(forest, ghost, quadrant, info.quadid[] + 0x1, info.treeid[] + 0x1, userdata)

        return
    end
    @cfunction($c, Cvoid, (Ptr{p8est_iter_volume_info_t}, Ptr{Cvoid}))
end

"""
    iterateforest(forest; kw...)

Execute the callbacks passed as keyword arguments for every volume, face,
edge, and corner of the rank-local forest.

The keyword arguments (`kw...`) for the iteration are:

  - `ghost = nothing`: the ghost layer associated for the mesh.  Used in
    `face`, `edge`, and `corner` callbacks for neighboring elements not
    rank local.
  - `volume = nothing`: Callback used for every volume (aka quadrant) of
    the local forest with the prototype
    `volume(forest, ghost, quadrant, quadid, treeid, userdata)`.
  - `face = nothing`: Not implemented yet.
  - `edge = nothing`: Not implemented yet.
  - `corner = nothing`: Not implemented yet.
  - `userdata = nothing`: User data passed to the callbacks.

See `@doc P4estTypes.P4est.p4est_iterate` and
`@doc P4estTypes.P4est.p8est_iterate` for more information about
the iteration.
"""
function iterateforest(
    forest::Pxest{4};
    ghost = nothing,
    volume = nothing,
    face = nothing,
    corner = nothing,
    userdata = nothing,
)
    if cfunction_closure && !isnothing(volume)
        volume_ = generate_volume_closure(Val(4), volume, forest, ghost, userdata)
    else
        volume_ = isnothing(volume) ? C_NULL : generate_volume_callback(Val(4))
    end

    ghost_ = isnothing(ghost) ? C_NULL : ghost
    face_::Ptr{Cvoid} = isnothing(face) ? C_NULL : error("Face iteration not implemented")
    corner_::Ptr{Cvoid} =
        isnothing(corner) ? C_NULL : error("Corner iteration not implemented")

    if cfunction_closure
        p4est_iterate(forest, ghost_, C_NULL, volume_, face_, corner_)
    else
        data = Ref((; forest, ghost, volume, face, corner, userdata))
        GC.@preserve data begin
            PointerWrapper(forest.pointer).user_pointer = pointer_from_objref(data)
            p4est_iterate(forest, ghost_, C_NULL, volume_, face_, corner_)
            PointerWrapper(forest.pointer).user_pointer = C_NULL
        end
    end

    return
end

function iterateforest(
    forest::Pxest{8};
    ghost = nothing,
    volume = nothing,
    face = nothing,
    edge = nothing,
    corner = nothing,
    userdata = nothing,
)
    if cfunction_closure && !isnothing(volume)
        volume_ = generate_volume_closure(Val(8), volume, forest, ghost, userdata)
    else
        volume_ = isnothing(volume) ? C_NULL : generate_volume_callback(Val(8))
    end

    ghost_ = isnothing(ghost) ? C_NULL : ghost
    face_::Ptr{Cvoid} = isnothing(face) ? C_NULL : error("Face iteration not implemented")
    edge_::Ptr{Cvoid} = isnothing(edge) ? C_NULL : error("Edge iteration not implemented")
    corner_::Ptr{Cvoid} =
        isnothing(corner) ? C_NULL : error("Corner iteration not implemented")

    if cfunction_closure
        p8est_iterate(forest, ghost_, C_NULL, volume_, face_, edge_, corner_)
    else
        data = Ref((; forest, ghost, volume, face, edge, corner, userdata))
        GC.@preserve data begin
            PointerWrapper(forest.pointer).user_pointer = pointer_from_objref(data)
            p8est_iterate(forest, ghost_, C_NULL, volume_, face_, edge_, corner_)
            PointerWrapper(forest.pointer).user_pointer = C_NULL
        end
    end

    return
end

function init_callback(forest, treeid, quadrant)
    data = unsafe_pointer_to_objref(pointer(PointerWrapper(forest).user_pointer))[]

    X = quadrantstyle(data.forest)
    Q = pxest_quadrant_t(Val(X))

    quadrant = QuadrantWrapper{X,Ptr{Q}}(quadrant)

    data.init(data.forest, treeid + 0x1, quadrant)

    return
end

@generated function generate_init_callback(::Val{X}) where {X}
    P = pxest_t(Val(X))
    Q = pxest_quadrant_t(Val(X))

    quote
        @cfunction(init_callback, Cvoid, (Ptr{$P}, p4est_topidx_t, Ptr{$Q}))
    end
end

function generate_init_closure(::Val{4}, init, forest)
    function c(::Ptr{p4est_t}, tid::p4est_topidx_t, q::Ptr{p4est_quadrant_t})::Cvoid
        quadrant = QuadrantWrapper{4,Ptr{p4est_quadrant_t}}(q)
        init(forest, tid + 0x1, quadrant)
        return
    end
    @cfunction($c, Cvoid, (Ptr{p4est_t}, p4est_topidx_t, Ptr{p4est_quadrant_t}))
end

function generate_init_closure(::Val{8}, init, forest)
    function c(::Ptr{p8est_t}, tid::p4est_topidx_t, q::Ptr{p8est_quadrant_t})::Cvoid
        quadrant = QuadrantWrapper{8,Ptr{p8est_quadrant_t}}(q)
        init(forest, tid + 0x1, quadrant)
        return
    end
    @cfunction($c, Cvoid, (Ptr{p8est_t}, p4est_topidx_t, Ptr{p8est_quadrant_t}))
end

function replace_callback(forest, treeid, num_outgoing, outgoing, num_incoming, incoming)
    data = unsafe_pointer_to_objref(pointer(PointerWrapper(forest).user_pointer))[]

    X = quadrantstyle(data.forest)
    Q = pxest_quadrant_t(Val(X))

    outgoing = unsafe_wrap(Array, outgoing, num_outgoing)
    outgoing = ntuple(i -> QuadrantWrapper{X,Ptr{Q}}(outgoing[i]), num_outgoing)

    incoming = unsafe_wrap(Array, incoming, num_incoming)
    incoming = ntuple(i -> QuadrantWrapper{X,Ptr{Q}}(incoming[i]), num_incoming)

    data.replace(data.forest, treeid + 0x1, outgoing, incoming)

    return
end

@generated function generate_replace_callback(::Val{X}) where {X}
    P = pxest_t(Val(X))
    Q = pxest_quadrant_t(Val(X))

    quote
        @cfunction(
            replace_callback,
            Cvoid,
            (Ptr{$P}, p4est_topidx_t, Cint, Ptr{Ptr{$Q}}, Cint, Ptr{Ptr{$Q}})
        )
    end
end

function generate_replace_closure(::Val{4}, replace, forest)
    function c(
        ::Ptr{p4est_t},
        tid::p4est_topidx_t,
        nout::Cint,
        out::Ptr{Ptr{p4est_quadrant_t}},
        ninx::Cint,
        inx::Ptr{Ptr{p4est_quadrant_t}},
    )::Cvoid
        outgoing =
            ntuple(i -> QuadrantWrapper{4,Ptr{p4est_quadrant_t}}(unsafe_load(out, i)), nout)
        incoming =
            ntuple(i -> QuadrantWrapper{4,Ptr{p4est_quadrant_t}}(unsafe_load(inx, i)), ninx)
        replace(forest, tid + 0x1, outgoing, incoming)
        return
    end
    @cfunction(
        $c,
        Cvoid,
        (
            Ptr{p4est_t},
            p4est_topidx_t,
            Cint,
            Ptr{Ptr{p4est_quadrant_t}},
            Cint,
            Ptr{Ptr{p4est_quadrant_t}},
        )
    )
end

function generate_replace_closure(::Val{8}, replace, forest)
    function c(
        ::Ptr{p8est_t},
        tid::p4est_topidx_t,
        nout::Cint,
        out::Ptr{Ptr{p8est_quadrant_t}},
        ninx::Cint,
        inx::Ptr{Ptr{p8est_quadrant_t}},
    )::Cvoid
        outgoing =
            ntuple(i -> QuadrantWrapper{8,Ptr{p8est_quadrant_t}}(unsafe_load(out, i)), nout)
        incoming =
            ntuple(i -> QuadrantWrapper{8,Ptr{p8est_quadrant_t}}(unsafe_load(inx, i)), ninx)
        replace(forest, tid + 0x1, outgoing, incoming)
        return
    end
    @cfunction(
        $c,
        Cvoid,
        (
            Ptr{p8est_t},
            p4est_topidx_t,
            Cint,
            Ptr{Ptr{p8est_quadrant_t}},
            Cint,
            Ptr{Ptr{p8est_quadrant_t}},
        )
    )
end

function coarsen_callback(forest, treeid, children)
    data = unsafe_pointer_to_objref(pointer(PointerWrapper(forest).user_pointer))[]

    X = quadrantstyle(data.forest)
    Q = pxest_quadrant_t(Val(X))

    children = unsafe_wrap(Array, children, X)
    children = ntuple(i -> QuadrantWrapper{X,Ptr{Q}}(children[i]), Val(X))
    return data.coarsen(data.forest, treeid + 0x1, children) ? one(Cint) : zero(Cint)
end

@generated function generate_coarsen_callback(::Val{X}) where {X}
    P = pxest_t(Val(X))
    Q = pxest_quadrant_t(Val(X))

    quote
        @cfunction(coarsen_callback, Cint, (Ptr{$P}, p4est_topidx_t, Ptr{Ptr{$Q}}))
    end
end

function generate_coarsen_closure(::Val{4}, coarsen, forest)
    function c(::Ptr{p4est_t}, tid::p4est_topidx_t, s::Ptr{Ptr{p4est_quadrant_t}})::Cint
        siblings =
            ntuple(i -> QuadrantWrapper{4,Ptr{p4est_quadrant_t}}(unsafe_load(s, i)), Val(4))
        return coarsen(forest, tid + 0x1, siblings)
    end
    @cfunction($c, Cint, (Ptr{p4est_t}, p4est_topidx_t, Ptr{Ptr{p4est_quadrant_t}}))
end

function generate_coarsen_closure(::Val{8}, coarsen, forest)
    function c(::Ptr{p8est_t}, tid::p4est_topidx_t, s::Ptr{Ptr{p8est_quadrant_t}})::Cint
        siblings =
            ntuple(i -> QuadrantWrapper{8,Ptr{p8est_quadrant_t}}(unsafe_load(s, i)), Val(8))
        return coarsen(forest, tid + 0x1, siblings)
    end
    @cfunction($c, Cint, (Ptr{p8est_t}, p4est_topidx_t, Ptr{Ptr{p8est_quadrant_t}}))
end

"""
    coarsen!(forest; coarsen = (_...) -> false, kw...)

Coarsen the quadrants of the forest determined by the `coarsen` callback.
The `coarsen(forest, treeid, siblings)` callback is called for each set
of sibling quadrants local to the rank that are eligible for coarsening.
If the callback returns `true` the `siblings` will coarsen into one quadrant
otherwise they will be untouched.

The other keyword arguments (`kw...`) for the coarsening are:

  - `recursive = false`: if `true` coarsening will be recursive otherwise each
    set of rank-local siblings will only be visited once.
  - `init = nothing`: callback function with prototype
    `init(forest, treeid, quadrant)` called for each quadrant created to
    initialized the user data.
  - `replace = nothing`: callback function with prototype
    `replace(forest, treeid, outgoing, incoming)` called for each set of
    `outgoing` quadrants with their associated `incoming` quadrant.  Note
    both `outgoing` and `incoming` are arrays with `eltype`
    [`QuadrantWrapper`](@ref).

See `@doc P4estTypes.P4est.p4est_coarsen_ext` and
`@doc P4estTypes.P4est.p8est_coarsen_ext` for more information about
the underlying p4est coarsening functions.
"""
function coarsen!(
    forest::Pxest{X};
    recursive = false,
    coarsen = (_...) -> false,
    init = nothing,
    replace = nothing,
) where {X}
    # Right now we do not support loading an orphan
    callback_orphans = false

    if cfunction_closure
        coarsen_ =
            isnothing(coarsen) ? C_NULL : generate_coarsen_closure(Val(X), coarsen, forest)
        init_ = isnothing(init) ? C_NULL : generate_init_closure(Val(X), init, forest)
        replace_ =
            isnothing(replace) ? C_NULL : generate_replace_closure(Val(X), replace, forest)

        (pxest_coarsen_ext(Val(X)))(
            forest,
            recursive,
            callback_orphans,
            coarsen_,
            init_,
            replace_,
        )
    else
        coarsen_ = isnothing(coarsen) ? C_NULL : generate_coarsen_callback(Val(X))
        init_ = isnothing(init) ? C_NULL : generate_init_callback(Val(X))
        replace_ = isnothing(replace) ? C_NULL : generate_replace_callback(Val(X))

        data = Ref((; forest, coarsen, init, replace))
        GC.@preserve data begin
            PointerWrapper(forest.pointer).user_pointer = pointer_from_objref(data)

            (pxest_coarsen_ext(Val(X)))(
                forest,
                recursive,
                callback_orphans,
                coarsen_,
                init_,
                replace_,
            )

            PointerWrapper(forest.pointer).user_pointer = C_NULL
        end
    end

    return
end

function refine_callback(forest, treeid, quadrant)
    data = unsafe_pointer_to_objref(pointer(PointerWrapper(forest).user_pointer))[]

    X = quadrantstyle(data.forest)
    Q = pxest_quadrant_t(Val(X))

    quadrant = QuadrantWrapper{X,Ptr{Q}}(quadrant)
    return data.refine(data.forest, treeid + 0x1, quadrant) ? one(Cint) : zero(Cint)
end

@generated function generate_refine_callback(::Val{X}) where {X}
    P = pxest_t(Val(X))
    Q = pxest_quadrant_t(Val(X))

    quote
        @cfunction(refine_callback, Cint, (Ptr{$P}, p4est_topidx_t, Ptr{$Q}))
    end
end

function generate_refine_closure(::Val{4}, refine, forest)
    function c(::Ptr{p4est_t}, tid::p4est_topidx_t, q::Ptr{p4est_quadrant_t})::Cint
        quadrant = QuadrantWrapper{4,Ptr{p4est_quadrant_t}}(q)
        return refine(forest, tid + 0x1, quadrant)
    end
    @cfunction($c, Cint, (Ptr{p4est_t}, p4est_topidx_t, Ptr{p4est_quadrant_t}))
end

function generate_refine_closure(::Val{8}, refine, forest)
    function c(::Ptr{p8est_t}, tid::p4est_topidx_t, q::Ptr{p8est_quadrant_t})::Cint
        quadrant = QuadrantWrapper{8,Ptr{p8est_quadrant_t}}(q)
        return refine(forest, tid + 0x1, quadrant)
    end
    @cfunction($c, Cint, (Ptr{p8est_t}, p4est_topidx_t, Ptr{p8est_quadrant_t}))
end

"""
    refine!(forest; refine = (_...) -> false, kw...)

Refine the quadrants of the forest determined by the `refine` callback.
The `refine(forest, treeid, quadrant)` callback is called for each quadrant
local to the rank. If the callback returns `true` the `quadrant` will
refine into multiple quadrants otherwise it will be untouched.

The other keyword arguments (`kw...`) for the refining are:

  - `recursive`: if `true` refining will be recursive otherwise each
    rank-local quadrant will only be visited once.
  - `maxlevel = -1`: the maximum level of refinement possible during this
    call.
  - `init = nothing`: callback function with prototype
    `init(forest, treeid, quadrant)` called for each quadrant created to
    initialized the user data.
  - `replace = nothing`: callback function with prototype
    `replace(forest, treeid, outgoing, incoming)` called for each
    `outgoing` quadrant with their associated `incoming` quadrants. Note both
    `outgoing` and `incoming` are arrays with `eltype` [`QuadrantWrapper`](@ref).

See `@doc P4estTypes.P4est.p4est_refine_ext` and
`@doc P4estTypes.P4est.p8est_refine_ext` for more information about
the underlying p4est refinement functions.
"""
function refine!(
    forest::Pxest{X};
    recursive = false,
    maxlevel = -1,
    refine = (_...) -> false,
    init = nothing,
    replace = nothing,
) where {X}
    if cfunction_closure
        refine_ =
            isnothing(refine) ? C_NULL : generate_refine_closure(Val(X), refine, forest)
        init_ = isnothing(init) ? C_NULL : generate_init_closure(Val(X), init, forest)
        replace_ =
            isnothing(replace) ? C_NULL : generate_replace_closure(Val(X), replace, forest)

        (pxest_refine_ext(Val(X)))(forest, recursive, maxlevel, refine_, init_, replace_)
    else
        refine_ = isnothing(refine) ? C_NULL : generate_refine_callback(Val(X))
        init_ = isnothing(init) ? C_NULL : generate_init_callback(Val(X))
        replace_ = isnothing(replace) ? C_NULL : generate_replace_callback(Val(X))

        data = Ref((; forest, refine, init, replace))

        GC.@preserve data begin
            PointerWrapper(forest.pointer).user_pointer = pointer_from_objref(data)

            (pxest_refine_ext(Val(X)))(
                forest,
                recursive,
                maxlevel,
                refine_,
                init_,
                replace_,
            )

            PointerWrapper(forest.pointer).user_pointer = C_NULL
        end
    end

    return
end

"""
    balance!(forest; kw...)

Enforce the two-to-one quadrant size constraint across the forest.  By default,
this constraint is enforced across faces, edges, and corners.

The keyword arguments (`kw...`) for the balancing are:

  - `connect`: type of constraint enforced which can take the values:

      + `P4estTypes.CONNECT_FULL(Val(4))`: enforce across face, and corner.
      + `P4estTypes.CONNECT_FULL(Val(8))`: enforce across face, edge, and corner.
      + `P4estTypes.CONNECT_FACE(Val(4))`: enforce across face.
      + `P4estTypes.CONNECT_FACE(Val(8))`: enforce across face.
      + `P4estTypes.CONNECT_EDGE(Val(8))`: enforce across face and edge.
      + `P4estTypes.CONNECT_CORNER(Val(4))`: enforce across face and corner.
      + `P4estTypes.CONNECT_CORNER(Val(8))`: enforce across face, edge, and corner.

  - `init = nothing`: callback function with prototype
    `init(forest, treeid, quadrant)` called for each quadrant created to
    initialized the user data.
  - `replace = nothing`: callback function with prototype
    `replace(forest, treeid, outgoing, incoming)` called for each
    `outgoing` quadrant with their associated `incoming` quadrants. Note both
    `outgoing` and `incoming` are arrays with `eltype` [`QuadrantWrapper`](@ref).

See `@doc P4estTypes.P4est.p4est_balance_ext` and
`@doc P4estTypes.P4est.p8est_balance_ext` for more information about
the underlying p4est balance functions.
"""
function balance!(
    forest::Pxest{X};
    connect = CONNECT_FULL(Val(X)),
    init = nothing,
    replace = nothing,
) where {X}
    if cfunction_closure
        init_ = isnothing(init) ? C_NULL : generate_init_closure(Val(X), init, forest)
        replace_ =
            isnothing(replace) ? C_NULL : generate_replace_closure(Val(X), replace, forest)

        (pxest_balance_ext(Val(X)))(forest, connect, init_, replace_)
    else
        init_ = isnothing(init) ? C_NULL : generate_init_callback(Val(X))
        replace_ = isnothing(replace) ? C_NULL : generate_replace_callback(Val(X))

        data = Ref((; forest, init, replace))

        GC.@preserve data begin
            PointerWrapper(forest.pointer).user_pointer = pointer_from_objref(data)

            (pxest_balance_ext(Val(X)))(forest, connect, init_, replace_)

            PointerWrapper(forest.pointer).user_pointer = C_NULL
        end
    end

    return
end

function weight_callback(forest, treeid, quadrant)
    data = unsafe_pointer_to_objref(pointer(PointerWrapper(forest).user_pointer))[]

    X = quadrantstyle(data.forest)
    Q = pxest_quadrant_t(Val(X))

    quadrant = QuadrantWrapper{X,Ptr{Q}}(quadrant)
    return data.weight(data.forest, treeid + 0x1, quadrant)
end

@generated function generate_weight_callback(::Val{X}) where {X}
    P = pxest_t(Val(X))
    Q = pxest_quadrant_t(Val(X))

    quote
        @cfunction(weight_callback, Cint, (Ptr{$P}, p4est_topidx_t, Ptr{$Q}))
    end
end

function generate_weight_closure(::Val{4}, weight, forest)
    function c(::Ptr{p4est_t}, tid::p4est_topidx_t, q::Ptr{p4est_quadrant_t})::Cint
        quadrant = QuadrantWrapper{4,Ptr{p4est_quadrant_t}}(q)
        return weight(forest, tid + 0x1, quadrant)
    end
    @cfunction($c, Cint, (Ptr{p4est_t}, p4est_topidx_t, Ptr{p4est_quadrant_t}))
end

function generate_weight_closure(::Val{8}, weight, forest)
    function c(::Ptr{p8est_t}, tid::p4est_topidx_t, q::Ptr{p8est_quadrant_t})::Cint
        quadrant = QuadrantWrapper{8,Ptr{p8est_quadrant_t}}(q)
        return weight(forest, tid + 0x1, quadrant)
    end
    @cfunction($c, Cint, (Ptr{p8est_t}, p4est_topidx_t, Ptr{p8est_quadrant_t}))
end

"""
    partition!(forest; kw...)

Partition the quadrants of the forest.  By default this will partition
the quadrants equally across the ranks of the forest.

By default sibling elements are split among the ranks.  This means they
cannot be coarsened with `coarsen!` and can cause MPI dependent coarsening.
If `allow_for_coarsening==true` then this is avoided by keeping sibling
quadrants on the same rank.

A `weight(forest, treeid, quadrant)` callback may provided (which gives the
`Float64` weight of each quadrant) for a weighted partitioning of the forest.

Alternatively, the forest may be partitioned to equally distribute the globally
numbered nodes via [`LNodes`](@ref).  This is done by setting `lnodes_degree`
to the node degree.  This requires the [`GhostLayer`](@ref) which if not
passed in `ghost` will be created.

The keyword arguments (`kw...`) for the partitioning are:

  - `ghost = nothing`: [`GhostLayer`](@ref) used when partitioning by
    [`LNodes`](@ref).
  - `lnodes_degree = nothing`: partition based on [`LNodes`](@ref) if this is
    set to the degree.
  - `allow_for_coarsening = false`: if `true` sibling groups that may be
    coarsened will be collect on the same rank.
  - `weight = nothing`: callback that give the `Float64` weight of each quadrant
    to perform a weighted partitioning.

See `@doc P4estTypes.P4est.p4est_partition_ext`,
`@doc P4estTypes.P4est.p8est_partition_ext`,
`@doc P4estTypes.P4est.p4est_partition_lnodes`,
and `@doc P4estTypes.P4est.p8est_partition_lnodes`, for more information about
the underlying p4est partition functions.
"""
function partition!(
    forest::Pxest{X};
    ghost = nothing,
    lnodes_degree = nothing,
    allow_for_coarsening = false,
    weight = nothing,
) where {X}
    if !isnothing(lnodes_degree)
        if isnothing(ghost)
            ghost = ghostlayer(forest)
        end
        (pxest_partition_lnodes(Val(X)))(forest, ghost, lnodes_degree, allow_for_coarsening)
    else
        if cfunction_closure
            weight_ =
                isnothing(weight) ? C_NULL : generate_weight_closure(Val(X), weight, forest)

            (pxest_partition_ext(Val(X)))(forest, allow_for_coarsening, weight_)
        else
            weight_ = isnothing(weight) ? C_NULL : generate_weight_callback(Val(X))

            data = Ref((; forest, weight))

            GC.@preserve data begin
                PointerWrapper(forest.pointer).user_pointer = pointer_from_objref(data)

                (pxest_partition_ext(Val(X)))(forest, allow_for_coarsening, weight_)

                PointerWrapper(forest.pointer).user_pointer = C_NULL
            end
        end
    end

    return
end

import Base: show
function Base.show(io::IO, forest::P4estTypes.Pxest{X}) where {X}
    print(io, "Forest{$X} with $(length(forest)) trees.")
end

function Base.show(io::IO, ::MIME{Symbol("text/plain")}, forest::P4estTypes.Pxest)
    print_tree(io, forest)
end

function Base.show(io::IO, tree::Tree{X}) where {X}
    print(io, "Tree{$X} with $(length(tree)) quadrants.")
end

function Base.show(io::IO, q::QuadrantWrapper{X}) where {X}
    print(io, "QuadrantWrapper{$X}: level $(level(q)), coordinates $(coordinates(q)).")
end

"""
    savevtk(prefix, forest; kw...)

Save the distributed forest-of-octrees (or quadtrees) `forest` to a set of
VTK files.

A `.vtu` file with the file name `prefix` is created per rank storing the
rank-local quadrants.  Further, `.pvtu` and `.visit` collection files
are created for ease of importing the mesh into Paraview and Visit,
respectively.

The keyword arguments (`kw...`) are:

  - `scale = 1.0`: a `scale < 1.0` places a visual gap between adjacent
    quadrants.
  - `writetree = true`: if `true` include the _zero-based_ tree id in VTK cell
    data.
  - `writelevel = true`: if `true` include the quadrant level in VTK cell data.
  - `writerank = true`: if `true` include the MPI rank in VTK cell data.
  - `wraprank = 0`: if `wraprank > 0` the MPI rank is stored modulo `wraprank`.
"""
function savevtk(
    prefix,
    forest::Pxest{X};
    scale = 1.0,
    writetree = true,
    writelevel = true,
    writerank = true,
    wraprank = false,
) where {X}
    context = (pxest_vtk_context_new(Val(X)))(forest, prefix)
    (pxest_vtk_context_set_geom(Val(X)))(context, C_NULL)
    (pxest_vtk_context_set_scale(Val(X)))(context, scale)
    (pxest_vtk_context_set_continuous(Val(X)))(context, true)
    context = (pxest_vtk_write_header(Val(X)))(context)
    if context == C_NULL
        error("pxest_vtk: Error writing header")
    end

    context = (pxest_vtk_write_cell_data(Val(X)))(
        context,
        writetree,
        writelevel,
        writerank,
        wraprank,
        0,
        0,
        C_NULL,
        C_NULL,
    )
    if context == C_NULL
        error("pxest_vtk: Error writing cell data")
    end

    if (pxest_vtk_write_footer(Val(X)))(context) != 0
        error("pxest_vtk: Error writing footer")
    end
end
