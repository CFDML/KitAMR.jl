"""
    LNodes{X,P}

Stores a parallel numbering of Lobatto nodes for a `Pxest{X}`.

# Fields
$(DocStringExtensions.FIELDS)

# See also
- [`lnodes`](@ref): a function used to construct `LNodes`
"""
mutable struct LNodes{X,P}
    """The pointer (of type `P`) can be a pointer to either a
    `P4estTypes.P4est.p4est_lnodes_t` or a
    `P4estTypes.P4est.p8est_lnodes_t`.  See the help
    documentation for these types for more information about the
    underlying p4est structures. """
    pointer::P
    """The MPI Communicator that includes the ranks participating in the
    lnods."""
    comm::MPI.Comm
    function LNodes{4}(pointer::Ptr{P4est.LibP4est.p4est_lnodes}, comm::MPI.Comm)
        nodes = new{4,typeof(pointer)}(pointer, comm)
        finalizer(nodes) do p
            p4est_lnodes_destroy(p.pointer)
            p.pointer = C_NULL
            return
        end
    end
    function LNodes{8}(pointer::Ptr{P4est.LibP4est.p8est_lnodes}, comm::MPI.Comm)
        nodes = new{8,typeof(pointer)}(pointer, comm)
        finalizer(nodes) do p
            p8est_lnodes_destroy(p.pointer)
            p.pointer = C_NULL
            return
        end
    end
end

"""
    lnodes(forest::Pxest; ghost = nothing, degree = 1)

Construct a parallel node numbering for the `forest`. If the
ghost layer `ghost` is not provided it will be constructed.

A `degree > 0` indicates that degree `degree` Lobotto nodes will be
constructed.  A `degree < 0` indicates that the boundary objects
(faces, edges, and corners) will be numbered.

See `@doc P4estTypes.P4est.p4est_lnodes_t` and
`@doc P4estTypes.P4est.p8est_lnodes_t` for a more detailed discussion
of the numbering based on `degree`.
"""
function lnodes(forest::Pxest{X}; ghost = nothing, degree = 1) where {X}
    if isnothing(ghost)
        ghost = ghostlayer(forest)
    end
    return LNodes{X}((pxest_lnodes_new(Val(X)))(forest, ghost, degree), comm(forest))
end

function Base.unsafe_convert(::Type{Ptr{p4est_lnodes}}, p::LNodes{4,Ptr{p4est_lnodes}})
    return p.pointer
end

function Base.unsafe_convert(::Type{Ptr{p8est_lnodes}}, p::LNodes{8,Ptr{p8est_lnodes}})
    return p.pointer
end

"""
    globalid(nodes::LNodes, localid::p4est_locidx_t)

Returns the global id associated with the node in `LNodes` with local id
`localid`.
"""
Base.@propagate_inbounds function globalid(nodes::LNodes, localid::p4est_locidx_t)
    pn = PointerWrapper(nodes.pointer)
    @boundscheck 0 ≤ localid < pn.num_local_nodes[]
    return localid < pn.owned_count[] ? pn.global_offset[] + localid :
           pn.nonlocal_nodes[localid-pn.owned_count[]+1]
end

"""
    unsafe_global_owned_count(nodes::LNodes)

Return an array containing the number of independent nodes owned by each
rank.

See `@doc P4estTypes.P4est.p4est_lnodes_t` and
`@doc P4estTypes.P4est.p8est_lnodes_t` for a more details.

Note, this unsafely wraps a C array.  So, you must ensure that the `nodes`
structure is preserved while using the return value.
"""
function unsafe_global_owned_count(nodes::LNodes)
    ns = PointerWrapper(nodes.pointer)
    mpisize = MPI.Comm_size(nodes.comm)

    global_owned_count = unsafe_wrap(
        Vector{p4est_locidx_t},
        pointer(ns.global_owned_count),
        (mpisize,),
        own = false,
    )

    return global_owned_count
end

"""
    unsafe_element_nodes(nodes::LNodes)

Return an array containing the unique continuous node number for each
local degree-of-freedom.

See `@doc P4estTypes.P4est.p4est_lnodes_t` and
`@doc P4estTypes.P4est.p8est_lnodes_t` for a more details.

Note, this unsafely wraps a C array.  So, you must ensure that the `nodes`
structure is preserved while using the return value.
"""
function unsafe_element_nodes(nodes::LNodes{X}) where {X}
    ns = PointerWrapper(nodes.pointer)

    element_nodes = unsafe_wrap(
        Matrix{p4est_locidx_t},
        pointer(ns.element_nodes),
        (ns.vnodes[], ns.num_local_elements[]),
        own = false,
    )

    if X == 4
        D = 2
    elseif X == 8
        D = 3
    else
        error("Not implemented")
    end

    return reshape(element_nodes, (ntuple(_ -> ns.degree[] + 1, D)..., :))
end

"""
    unsafe_face_code(nodes::LNodes)

Return an array containing the face code for each quadrant of the mesh. The
face code indicates which faces and edges of the quadrant are hanging.

See the p4est functions `p4est_lnodes_decode` and `p8est_lnodes_decode` to
determine how to decode the face code.

Note, this unsafely wraps a C array.  So, you must ensure that the `nodes`
structure is preserved while using the return value.
"""
function unsafe_face_code(nodes::LNodes)
    ns = PointerWrapper(nodes.pointer)

    ptr = pointer(ns.face_code)

    element_nodes =
        unsafe_wrap(Vector{eltype(ptr)}, ptr, (ns.num_local_elements[],), own = false)

    return element_nodes
end

"""
    sharers(nodes::LNodes)

Returns a `Dict` mapping the neighboring rank with the global ids of the
nodes shared between it and the local rank.

Note, this `Dict` does not include an entry for the local rank.
"""
function sharers(nodes::LNodes)
    d = Dict{Cint,Set{p4est_gloidx_t}}()
    rank = MPI.Comm_rank(nodes.comm)

    pn = PointerWrapper(nodes.pointer)
    sa = PointerWrapper(Ptr{p4est_lnodes_rank}(pointer(pn.sharers.array)))
    for n = 1:pn.sharers.elem_count[]
        if sa[n].rank == rank
            continue
        end

        s = get!(d, sa[n].rank[]) do
            Set{p4est_gloidx_t}()
        end

        for m = 1:sa[n].shared_nodes.elem_count
            lid = unsafe_load(Ptr{p4est_locidx_t}(sa[n].shared_nodes.array), m)
            gid = globalid(nodes, lid)
            push!(s, gid)
        end
    end

    return d
end
