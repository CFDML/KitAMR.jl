# ---------------------------------------------------------------------------
# Closure-free C callbacks (Apple Silicon / aarch64 does not support closure
# `@cfunction`). A fixed top-level trampoline is wrapped once with `@cfunction`
# and the runtime-chosen Julia callable is recovered from a pointer the C side
# can reach:
#   * weight callbacks have no user_data slot, so the callable is threaded through
#     the forest `user_pointer` inside an `AMR_weight_context` (which also keeps
#     the object `user_pointer` originally referenced, e.g. `kinfo`/`ka`). The
#     trampoline restores that original object around each call, and the partition
#     wrappers restore `user_pointer` once partitioning finishes. This is reentrant
#     — the state lives on the forest, not in a module-global Ref.
#   * init callbacks (forest creation) still use a module-level Ref read by a fixed
#     trampoline; the previous Ref value is restored after the call.
#   * volume/face/corner iterate callbacks carry the callable (plus data_type
#     and the caller's user_data) in `AMR_iterate_context`, passed through the
#     p4est_iterate `user_data` channel and recovered in a top-level trampoline.
#     The forest `user_pointer` is never disturbed.
# ---------------------------------------------------------------------------

"""
Context threaded through the forest `user_pointer` so a single fixed trampoline
`@cfunction` can dispatch to a runtime-chosen (possibly closure) init function at
forest creation, without a module-global `Ref`. Mirrors [`AMR_weight_context`]:
`user_pointer` keeps the object originally passed as the forest `user_pointer`
(e.g. `kinfo`), which the trampoline restores around each init call; the
`AMR_4est_new` wrappers restore `user_pointer` once creation finishes.
"""
mutable struct AMR_init_context
    f::Function
    user_pointer::Ptr{Cvoid}
end

"""
Context threaded through the forest `user_pointer` so a single fixed trampoline
`@cfunction` can dispatch to a runtime-chosen (possibly closure) weight function
without a module-global `Ref` (which was non-reentrant: only one partition could be
in flight at a time, and nesting would clobber it). `user_pointer` keeps whatever
the forest `user_pointer` originally referenced (e.g. `kinfo`/`ka`); the trampoline
temporarily restores it around each weight call so weight functions that read
`user_pointer` keep working.
"""
mutable struct AMR_weight_context
    f::Function
    user_pointer::Ptr{Cvoid}
end

function _amr_weight_cb_p4est(p4est::Ptr{p4est_t}, which_tree::p4est_topidx_t, quadrant::Ptr{p4est_quadrant_t})::Cint
    fp = PointerWrapper(p4est)
    ctx = unsafe_pointer_to_objref(pointer(fp.user_pointer))::AMR_weight_context
    GC.@preserve ctx begin
        fp.user_pointer = ctx.user_pointer          # expose the original object to the weight fn
        w = Cint(ctx.f(p4est, which_tree, quadrant))
        fp.user_pointer = pointer_from_objref(ctx)  # restore ctx for the next quadrant
    end
    return w
end
function _amr_weight_cb_p8est(p4est::Ptr{p8est_t}, which_tree::p4est_topidx_t, quadrant::Ptr{p8est_quadrant_t})::Cint
    fp = PointerWrapper(p4est)
    ctx = unsafe_pointer_to_objref(pointer(fp.user_pointer))::AMR_weight_context
    GC.@preserve ctx begin
        fp.user_pointer = ctx.user_pointer
        w = Cint(ctx.f(p4est, which_tree, quadrant))
        fp.user_pointer = pointer_from_objref(ctx)
    end
    return w
end
function _amr_init_cb_p4est(p4est::Ptr{p4est_t}, which_tree::p4est_topidx_t, quadrant::Ptr{p4est_quadrant_t})::Cvoid
    fp = PointerWrapper(p4est)
    ctx = unsafe_pointer_to_objref(pointer(fp.user_pointer))::AMR_init_context
    GC.@preserve ctx begin
        fp.user_pointer = ctx.user_pointer          # expose the original object to the init fn
        ctx.f(p4est, which_tree, quadrant)
        fp.user_pointer = pointer_from_objref(ctx)  # restore ctx for the next quadrant
    end
    return nothing
end
function _amr_init_cb_p8est(p4est::Ptr{p8est_t}, which_tree::p4est_topidx_t, quadrant::Ptr{p8est_quadrant_t})::Cvoid
    fp = PointerWrapper(p4est)
    ctx = unsafe_pointer_to_objref(pointer(fp.user_pointer))::AMR_init_context
    GC.@preserve ctx begin
        fp.user_pointer = ctx.user_pointer
        ctx.f(p4est, which_tree, quadrant)
        fp.user_pointer = pointer_from_objref(ctx)
    end
    return nothing
end

"""
Bundles a Julia callable together with the per-call `user_data` (and, for volume
iteration, the quadrant `data_type`) so it can be passed into a p4est_iterate C
callback through the `user_data` channel without forming a closure.
"""
mutable struct AMR_iterate_context
    f::Function
    data_type::DataType
    user_data::Ptr{Cvoid}
end

function _amr_volume_iter_cb_p4est(info::Ptr{p4est_iter_volume_info}, data::Ptr{Nothing})::Cvoid
    GC.@preserve info begin
        ctx = unsafe_pointer_to_objref(data)::AMR_iterate_context
        ip = PointerWrapper(info)
        dp = PointerWrapper(ctx.data_type, ip.quad.p.user_data[])
        GC.@preserve ip dp ctx.f(ip, ctx.user_data, dp)
    end
    return nothing
end
function _amr_volume_iter_cb_p8est(info::Ptr{p8est_iter_volume_info}, data::Ptr{Nothing})::Cvoid
    GC.@preserve info begin
        ctx = unsafe_pointer_to_objref(data)::AMR_iterate_context
        ip = PointerWrapper(info)
        dp = PointerWrapper(ctx.data_type, ip.quad.p.user_data[])
        GC.@preserve ip dp ctx.f(ip, ctx.user_data, dp)
    end
    return nothing
end
function _amr_corner_iter_cb_p4est(info::Ptr{p4est_iter_corner_info}, data::Ptr{Nothing})::Cvoid
    GC.@preserve info begin
        ctx = unsafe_pointer_to_objref(data)::AMR_iterate_context
        ip = PointerWrapper(info)
        GC.@preserve ip ctx.f(ip, ctx.user_data)
    end
    return nothing
end
function _amr_face_iter_cb_p4est(info::Ptr{p4est_iter_face_info}, data::Ptr{Nothing})::Cvoid
    GC.@preserve info begin
        ctx = unsafe_pointer_to_objref(data)::AMR_iterate_context
        ip = PointerWrapper(info)
        GC.@preserve ip ctx.f(ip, ctx.user_data)
    end
    return nothing
end
function _amr_face_iter_cb_p8est(info::Ptr{p8est_iter_face_info}, data::Ptr{Nothing})::Cvoid
    GC.@preserve info begin
        ctx = unsafe_pointer_to_objref(data)::AMR_iterate_context
        ip = PointerWrapper(info)
        GC.@preserve ip ctx.f(ip, ctx.user_data)
    end
    return nothing
end

function quad_to_cell(
    fp::PointerWrapper{p4est_t},
    treeid,
    qp::PointerWrapper{p4est_quadrant_t},
)
    quad = zeros(3)
    mid = zeros(3)
    GC.@preserve quad mid fp qp treeid begin
        halflength = Cint(round(P4EST_QUADRANT_LEN(qp.level[]) / 2))
        p4est_qcoord_to_vertex(pointer(fp.connectivity), treeid, qp.x[], qp.y[], quad)
        p4est_qcoord_to_vertex(
            pointer(fp.connectivity),
            treeid,
            qp.x[] + halflength,
            qp.y[] + halflength,
            mid,
        )
        ds = 2 * (mid - quad)
    end
    (ds[1:2], mid[1:2])
end
function quad_to_cell(
    fp::PointerWrapper{p8est_t},
    treeid,
    qp::PointerWrapper{p8est_quadrant_t},
)
    quad = zeros(3)
    mid = zeros(3)
    GC.@preserve quad mid fp qp treeid begin
        halflength = Cint(round(P8EST_QUADRANT_LEN(qp.level[]) / 2))
        p8est_qcoord_to_vertex(
            pointer(fp.connectivity),
            treeid,
            qp.x[],
            qp.y[],
            qp.z[],
            quad,
        )
        p8est_qcoord_to_vertex(
            pointer(fp.connectivity),
            treeid,
            qp.x[] + halflength,
            qp.y[] + halflength,
            qp.z[] + halflength,
            mid,
        )
        ds = 2. * (mid - quad)
    end
    (ds, mid)
end

function iPointerWrapper(p::Ptr{T},i::Integer) where{T}
    return PointerWrapper(p + sizeof(T) * i)
end
function iPointerWrapper(pw::PointerWrapper{T},i::Integer) where{T}
    return PointerWrapper(pointer(pw)+i*sizeof(T))
end
function iPointerWrapper(p::Ptr{sc_array_t},::Type{T},i::Integer) where{T}
    return PointerWrapper(T,pointer(PointerWrapper(p).array)+i*sizeof(T))
end
function iPointerWrapper(p::Ptr{sc_array_t},::Type{Ptr{T}},i::Integer) where{T} # to avoid the dereference of the Pointer-type
    return PointerWrapper(Ptr{Ptr{T}}(pointer(PointerWrapper(p).array)+i*sizeof(Ptr{T})))
end
function iPointerWrapper(pw::PointerWrapper{sc_array_t},::Type{T},i::Integer) where{T}
    return PointerWrapper(T,pointer(pw.array)+i*sizeof(T))
end
function iPointerWrapper(pw::PointerWrapper{sc_array_t},::Type{Ptr{T}},i::Integer) where{T}
    return PointerWrapper(Ptr{Ptr{T}}(pointer(pw.array)+i*sizeof(Ptr{T})))
end
function iPointerWrapper(pw::PointerWrapper{NTuple{N,T}},::Type{T},i::Integer) where{N,T}
    return PointerWrapper(T,pointer(pw)+i*sizeof(T))
end
function iPointerWrapper(pw::PointerWrapper{NTuple{N,Ptr{T}}},::Type{Ptr{T}},i::Integer) where{N,T}
    return PointerWrapper(Ptr{Ptr{T}}(pointer(pw)+i*sizeof(Ptr{T})))
end

function local_quadid(p4est::PointerWrapper{p4est_t},treeid::Integer,quadid::Integer)
    tp = iPointerWrapper(p4est.trees, p4est_tree_t, treeid)
    return tp.quadrants_offset[] + quadid
end
function local_quadid(ip::PointerWrapper{p4est_iter_volume_info_t})
    tp = iPointerWrapper(ip.p4est.trees, p4est_tree_t, ip.treeid[])
    return tp.quadrants_offset[] + ip.quadid[]
end
function local_quadid(ip::PointerWrapper{p4est_iter_corner_info_t},side::PointerWrapper{p4est_iter_corner_side_t})
    tp = iPointerWrapper(ip.p4est.trees, p4est_tree_t, side.treeid[])
    return tp.quadrants_offset[] + side.quadid[]
end
function local_quadid(ip::PointerWrapper{p8est_iter_volume_info_t})
    tp = iPointerWrapper(ip.p4est.trees, p8est_tree_t, ip.treeid[])
    return tp.quadrants_offset[] + ip.quadid[]
end
function local_quadid(ip::PointerWrapper{p8est_iter_corner_info_t},side::PointerWrapper{p8est_iter_corner_side_t})
    tp = iPointerWrapper(ip.p4est.trees, p4est_tree_t, side.treeid[])
    return tp.quadrants_offset[] + side.quadid[]
end
function global_quadid(ip::PW_pxest_iter_volume_info_t)
    gfq = unsafe_wrap(
            Vector{Int},
            pointer(ip.p4est.global_first_quadrant),
            MPI.Comm_size(MPI.COMM_WORLD) + 1,
        )
    return local_quadid(ip) + gfq[MPI.Comm_rank(MPI.COMM_WORLD)+1]
end


function AMR_4est_new(
    comm::MPI.Comm,
    conn::Ptr{p4est_connectivity},
    min_quads::Integer,
    min_level::Integer,
    uniform::Integer,
    ::Type{T},
    init_fn::Function,
    user_pointer::Ptr,
) where {T}
    orig = Ptr{Cvoid}(user_pointer)
    ctx = AMR_init_context(init_fn, orig)
    GC.@preserve comm conn init_fn ctx begin
        p4est = p4est_new_ext(
            comm,
            conn,
            min_quads,
            min_level,
            uniform,
            sizeof(T),
            @cfunction(_amr_init_cb_p4est, Cvoid, (Ptr{p4est_t}, p4est_topidx_t, Ptr{p4est_quadrant_t})),
            pointer_from_objref(ctx),
        )
        PointerWrapper(p4est).user_pointer = orig   # restore original user_pointer
        return p4est
    end
end
function AMR_4est_new(
    comm::MPI.Comm,
    conn::Ptr{p8est_connectivity},
    min_quads::Integer,
    min_level::Integer,
    uniform::Integer,
    ::Type{T},
    init_fn::Function,
    user_pointer::Ptr,
) where {T}
    orig = Ptr{Cvoid}(user_pointer)
    ctx = AMR_init_context(init_fn, orig)
    GC.@preserve comm conn init_fn ctx begin
        p8est = p8est_new_ext(
            comm,
            conn,
            min_quads,
            min_level,
            uniform,
            sizeof(T),
            @cfunction(_amr_init_cb_p8est, Cvoid, (Ptr{p8est_t}, p4est_topidx_t, Ptr{p8est_quadrant_t})),
            pointer_from_objref(ctx),
        )
        PointerWrapper(p8est).user_pointer = orig   # restore original user_pointer
        return p8est
    end
end

function AMR_4est_new(comm::MPI.Comm, conn::Ptr{p4est_connectivity}, ::Type{T}) where {T}
    GC.@preserve comm conn p4est_new_ext(comm, conn, 1, 0, 1, sizeof(T), C_NULL, C_NULL)
end
function AMR_4est_new(comm::MPI.Comm, conn::Ptr{p8est_connectivity}, ::Type{T}) where {T}
    GC.@preserve comm conn p8est_new_ext(comm, conn, 1, 0, 1, sizeof(T), C_NULL, C_NULL)
end

function AMR_4est_new(
    comm::MPI.Comm,
    conn::Ptr{p4est_connectivity},
    ::Type{T},
    init_fn::Function,
    user_pointer::Ptr,
) where {T}
    orig = Ptr{Cvoid}(user_pointer)
    ctx = AMR_init_context(init_fn, orig)
    GC.@preserve comm conn init_fn ctx begin
        p4est = p4est_new_ext(
            comm,
            conn,
            1,
            0,
            1,
            sizeof(T),
            @cfunction(_amr_init_cb_p4est, Cvoid, (Ptr{p4est_t}, p4est_topidx_t, Ptr{p4est_quadrant_t})),
            pointer_from_objref(ctx),
        )
        PointerWrapper(p4est).user_pointer = orig   # restore original user_pointer
        return p4est
    end
end
function AMR_4est_new(
    comm::MPI.Comm,
    conn::Ptr{p8est_connectivity},
    ::Type{T},
    init_fn::Function,
    user_pointer::Ptr,
) where {T}
    orig = Ptr{Cvoid}(user_pointer)
    ctx = AMR_init_context(init_fn, orig)
    GC.@preserve comm conn init_fn ctx begin
        p8est = p8est_new_ext(
            comm,
            conn,
            1,
            0,
            1,
            sizeof(T),
            @cfunction(_amr_init_cb_p8est, Cvoid, (Ptr{p8est_t}, p4est_topidx_t, Ptr{p8est_quadrant_t})),
            pointer_from_objref(ctx),
        )
        PointerWrapper(p8est).user_pointer = orig   # restore original user_pointer
        return p8est
    end
end

function AMR_4est_new(
    comm::MPI.Comm,
    conn::Ptr{p4est_connectivity},
    ::Type{T},
    user_pointer::Ptr,
) where {T}
    GC.@preserve comm conn user_pointer p4est_new_ext(
        comm,
        conn,
        1,
        0,
        1,
        sizeof(T),
        C_NULL,
        user_pointer,
    )
end
function AMR_4est_new(
    comm::MPI.Comm,
    conn::Ptr{p8est_connectivity},
    ::Type{T},
    user_pointer::Ptr,
) where {T}
    GC.@preserve comm conn user_pointer p8est_new_ext(
        comm,
        conn,
        1,
        0,
        1,
        sizeof(T),
        C_NULL,
        user_pointer,
    )
end

function AMR_4est_new(
    comm::MPI.Comm,
    conn::Ptr{p4est_connectivity},
    ::Type{T},
    init_fn::Function,
) where {T}
    orig = Ptr{Cvoid}(C_NULL)
    ctx = AMR_init_context(init_fn, orig)
    GC.@preserve comm conn init_fn ctx begin
        p4est = p4est_new_ext(
            comm,
            conn,
            1,
            0,
            1,
            sizeof(T),
            @cfunction(_amr_init_cb_p4est, Cvoid, (Ptr{p4est_t}, p4est_topidx_t, Ptr{p4est_quadrant_t})),
            pointer_from_objref(ctx),
        )
        PointerWrapper(p4est).user_pointer = orig   # restore original user_pointer (C_NULL)
        return p4est
    end
end
function AMR_4est_new(
    comm::MPI.Comm,
    conn::Ptr{p8est_connectivity},
    ::Type{T},
    init_fn::Function,
) where {T}
    orig = Ptr{Cvoid}(C_NULL)
    ctx = AMR_init_context(init_fn, orig)
    GC.@preserve comm conn init_fn ctx begin
        p8est = p8est_new_ext(
            comm,
            conn,
            1,
            0,
            1,
            sizeof(T),
            @cfunction(_amr_init_cb_p8est, Cvoid, (Ptr{p8est_t}, p4est_topidx_t, Ptr{p8est_quadrant_t})),
            pointer_from_objref(ctx),
        )
        PointerWrapper(p8est).user_pointer = orig   # restore original user_pointer (C_NULL)
        return p8est
    end
end

function AMR_face_iterate(
    info::P_pxest_iter_face_info_t,
    data::Ptr{Nothing},
    kernel::Function,
)
    GC.@preserve info data kernel begin
        ip = PointerWrapper(info)
        GC.@preserve ip kernel(ip, data)
    end
    return nothing
end

"""
$(TYPEDSIGNATURES)
"""
function AMR_ghost_new(p4est::Ptr{p4est_t})
    GC.@preserve p4est p4est_ghost_new(p4est, P4EST_CONNECT_FULL)
end
"""
$(TYPEDSIGNATURES)
"""
function AMR_ghost_new(p4est::Ptr{p8est_t})
    GC.@preserve p4est p8est_ghost_new(p4est, P8EST_CONNECT_FULL)
end

function AMR_partition(p4est::Ptr{p4est_t},weight_fn::T=C_NULL) where{T<:Union{Ptr{Nothing},Base.CFunction}}
    p4est_partition(p4est, 0, weight_fn)
end
function AMR_partition(weight_fn::Function,p4est::Ptr{p4est_t})
    fp = PointerWrapper(p4est)
    orig = Ptr{Cvoid}(pointer(fp.user_pointer))            # original user_pointer (e.g. kinfo/ka)
    ctx = AMR_weight_context(weight_fn, orig)
    GC.@preserve weight_fn ctx begin
        fp.user_pointer = pointer_from_objref(ctx)
        try
            AMR_partition(p4est,
                @cfunction(_amr_weight_cb_p4est,Cint,(Ptr{p4est_t},p4est_topidx_t,Ptr{p4est_quadrant_t})))
        finally
            fp.user_pointer = orig                         # restore original user_pointer
        end
    end
end
function AMR_partition(p4est::Ptr{p8est_t},weight_fn::T=C_NULL) where{T<:Union{Ptr{Nothing},Base.CFunction}}
    p8est_partition(p4est, 0, weight_fn)
end
function AMR_partition(weight_fn::Function,p4est::Ptr{p8est_t})
    fp = PointerWrapper(p4est)
    orig = Ptr{Cvoid}(pointer(fp.user_pointer))            # original user_pointer (e.g. kinfo/ka)
    ctx = AMR_weight_context(weight_fn, orig)
    GC.@preserve weight_fn ctx begin
        fp.user_pointer = pointer_from_objref(ctx)
        try
            AMR_partition(p4est,
                @cfunction(_amr_weight_cb_p8est,Cint,(Ptr{p8est_t},p4est_topidx_t,Ptr{p8est_quadrant_t})))
        finally
            fp.user_pointer = orig                         # restore original user_pointer
        end
    end
end

"""
$(TYPEDSIGNATURES)
"""
function AMR_mesh_new(p4est::Ptr{p4est_t},ghost::Ptr{p4est_ghost_t})
    GC.@preserve p4est ghost p4est_mesh_new_ext(p4est,ghost,1,1,P4EST_CONNECT_FULL)
end
"""
$(TYPEDSIGNATURES)
"""
function AMR_mesh_new(p4est::Ptr{p8est_t},ghost::Ptr{p8est_ghost_t})
    GC.@preserve p4est ghost p8est_mesh_new_ext(p4est,ghost,1,1,P8EST_CONNECT_FULL)
end


function unsafe_wrap_sc(::Type{T}, sc_array_ptr::Ptr{sc_array}) where {T}
    sc_array_obj = unsafe_load(sc_array_ptr)
    return unsafe_wrap_sc(T, sc_array_obj)
end

function unsafe_wrap_sc(::Type{T}, sc_array_obj::sc_array) where {T}
    elem_count = sc_array_obj.elem_count
    array = sc_array_obj.array
    return unsafe_wrap(Vector{T}, Ptr{T}(array), elem_count)
end
function unsafe_wrap_sc(::Type{T}, sc_array_pw::PointerWrapper{sc_array}) where {T}
    elem_count = sc_array_pw.elem_count[]
    array = sc_array_pw.array
    return unsafe_wrap(Vector{T}, Ptr{T}(pointer(array)), elem_count)
end

function AMR_volume_iterate(f::Function,forest::Ptr{p4est_t};ghost=C_NULL,user_data=C_NULL,data_type = P4estPsData)
    ctx = AMR_iterate_context(f, data_type, Ptr{Cvoid}(user_data))
    GC.@preserve f ctx forest ghost p4est_iterate(
        forest,
        ghost,
        pointer_from_objref(ctx),
        @cfunction(_amr_volume_iter_cb_p4est,Cvoid, (Ptr{p4est_iter_volume_info}, Ptr{Nothing})),
        C_NULL,
        C_NULL,
    )
end
function AMR_volume_iterate(f::Function,forest::Ptr{p8est_t};ghost=C_NULL,user_data=C_NULL,data_type = P4estPsData)
    ctx = AMR_iterate_context(f, data_type, Ptr{Cvoid}(user_data))
    GC.@preserve f ctx forest ghost p8est_iterate(
        forest,
        ghost,
        pointer_from_objref(ctx),
        @cfunction(_amr_volume_iter_cb_p8est,Cvoid, (Ptr{p8est_iter_volume_info}, Ptr{Nothing})),
        C_NULL,
        C_NULL,
        C_NULL,
    )
end
function AMR_corner_iterate(f::Function,forest::Ptr{p4est_t};ghost=C_NULL,user_data=C_NULL)
    ctx = AMR_iterate_context(f, Nothing, Ptr{Cvoid}(user_data))
    GC.@preserve f ctx forest ghost p4est_iterate(
        forest,
        ghost,
        pointer_from_objref(ctx),
        C_NULL,
        C_NULL,
        @cfunction(_amr_corner_iter_cb_p4est,Cvoid, (Ptr{p4est_iter_corner_info}, Ptr{Nothing})),
    )
end
function AMR_face_iterate(f::Function,forest::Ptr{p4est_t};ghost=C_NULL,user_data=C_NULL)
    ctx = AMR_iterate_context(f, Nothing, Ptr{Cvoid}(user_data))
    GC.@preserve f ctx forest ghost p4est_iterate(
        forest,
        ghost,
        pointer_from_objref(ctx),
        C_NULL,
        @cfunction(_amr_face_iter_cb_p4est,Cvoid, (Ptr{p4est_iter_face_info}, Ptr{Nothing})),
        C_NULL,
    )
end
function AMR_face_iterate(f::Function,forest::Ptr{p8est_t};ghost=C_NULL,user_data=C_NULL)
    ctx = AMR_iterate_context(f, Nothing, Ptr{Cvoid}(user_data))
    GC.@preserve f ctx forest ghost p8est_iterate( # volume, face, edge, corner
        forest,
        ghost,
        pointer_from_objref(ctx),
        C_NULL,
        @cfunction(_amr_face_iter_cb_p8est,Cvoid, (Ptr{p8est_iter_face_info}, Ptr{Nothing})),
        C_NULL,
        C_NULL,
    )
end