function get_index_ptr(p::PointerWrapper{T1}, ::Type{T2}, i::Integer) where {T1,T2}
    return Ptr{T2}(pointer(p) + i * sizeof(T2))
end
function get_index_ptr(
    array::Ptr{T},
    t::DataType,
    i::Union{Int,UInt,Cint},
) where {T<:sc_array}
    sc = PointerWrapper(array)
    return Ptr{t}(pointer(sc.array) + i * sc.elem_size[])
end
function get_index_ptr(
    pw_sc_array::PointerWrapper{sc_array_t},
    t::DataType,
    i::Union{Int,UInt,Cint},
)
    return Ptr{t}(pointer(pw_sc_array.array) + i * pw_sc_array.elem_size[])
end
function get_index_ptr(ghost_datas::Ptr{T}, t::DataType, i::Union{Int,UInt,Cint}) where {T}
    if t != T
        @warn "get_index_ptr: type mismatch"
    end
    return Ptr{t}(ghost_datas + i * sizeof(t))
end
function get_index_pw(ptr::Ptr{T}, i::Integer) where {T}
    return PointerWrapper(ptr + sizeof(T) * i)
end
function get_index_pw(array, t::DataType, i::Union{Int,UInt,Cint})
    return PointerWrapper(get_index_ptr(array, t, i))
end
function get_index_pw(p::PointerWrapper{T}, i::Integer) where {T}
    return PointerWrapper(pointer(p) + sizeof(T) * i)
end
function get_midpoint(forest::Ptr{p4est_t}, treeid, quadrant::Ptr{p4est_quadrant_t})
    fp = PointerWrapper(forest)
    qp = PointerWrapper(quadrant)
    get_midpoint(fp, treeid, qp)
end
function get_midpoint(
    fp::PointerWrapper{p4est_t},
    treeid,
    qp::PointerWrapper{p4est_quadrant_t},
)
    midpoint = zeros(3)
    GC.@preserve midpoint fp qp treeid begin
        halflength = P4EST_QUADRANT_LEN(qp.level[]) / 2
        x = qp.x[] + halflength
        y = qp.y[] + halflength
        p4est_qcoord_to_vertex(pointer(fp.connectivity), treeid, x, y, midpoint)
    end
    midpoint[1:2]
end
function get_midpoint_ds(forest::Ptr{p4est_t}, treeid, quadrant::Ptr{p4est_quadrant_t})
    fp = PointerWrapper(forest)
    qp = PointerWrapper(quadrant)
    get_midpoint_ds(fp, treeid, qp)
end
function get_midpoint_ds(
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
    (ds[1:DIM], mid[1:DIM])
end
function ghost_data_alloc(t::DataType, ghost::Ptr{p4est_ghost_t})
    return GC.@preserve t ghost Ptr{t}(
        sc_malloc(-1, sizeof(t) * PointerWrapper(ghost).ghosts.elem_count[]),
    )
end
function ghost_data_alloc(N::Int, ghost::Ptr{p4est_ghost_t})
    return GC.@preserve N ghost Ptr{Cdouble}(
        sc_malloc(-1, sizeof(Cdouble) * N * PointerWrapper(ghost).ghosts.elem_count[]),
    )
end
function ghost_data_pw(ghost_data::Ptr{T}, qid::Integer) where {T}
    return PointerWrapper(T, ghost_data_ptr(T, qid, ghost_data))
end
function ghost_data_ptr(t::DataType, qid::Integer, ghost_data::Ptr{T}) where {T}
    if T != t
        @warn "ghost_data_ptr: type mismatch"
    end
    return Ptr{t}(ghost_data + qid * sizeof(t))
end
function ghost_data_ptr(N::Int, ghost_data::Ptr{T}, qid::Integer) where {T}
    return ghost_data + qid * sizeof(T) * N
end
function get_neighbor_data_ptr(
    ip::PointerWrapper{p4est_iter_volume_info_t},
    data::Reconstruct,
    dir::Integer,
)
    quadid = quadid(ip)
    global_data = data.global_data
    ghost_datas = data.ghost_datas
    get_neighbor_data_ptr(
        pointer(ip.p4est),
        global_data.ghost,
        global_data.mesh,
        quadid,
        ghost_datas,
        eltype(ghost_datas),
        dir,
    )
end
function get_neighbor_data_ptr(
    p4est::Ptr{p4est_t},
    ghost::Ptr{p4est_ghost_t},
    mesh::Ptr{p4est_mesh_t},
    quadid::p4est_locidx_t,
    ghost_data::Ptr,
    ::Type{t},
    direction::Union{Int,UInt},
) where {t}  # 0: left, 1: right, 2: bottom, 3: top 
    neighbor_quads = sc_array_new(sizeof(Ptr{p4est_quadrant_t}))
    neighbor_encs = sc_array_new(sizeof(Cint))
    neighbor_qid = sc_array_new(sizeof(Cint))
    GC.@preserve p4est ghost mesh quadid direction neighbor_encs neighbor_qid neighbor_quads t begin
        p4est_mesh_get_neighbors(
            p4est,
            ghost,
            mesh,
            quadid,
            direction,
            neighbor_quads,
            neighbor_encs,
            neighbor_qid,
        )
        scpw = PointerWrapper(neighbor_qid)
        if scpw.elem_count[] > 0
            enc = PointerWrapper(get_index_ptr(neighbor_encs, Cint, 0))[]
            if enc < 0
                qid = PointerWrapper(get_index_ptr(scpw, Cint, 0))[]
                sc_array_destroy(neighbor_quads)
                sc_array_destroy(neighbor_encs)
                sc_array_destroy(neighbor_qid)
                return ghost_data_ptr(t, qid, ghost_data)
            else
                pq = PointerWrapper(
                    PointerWrapper(
                        get_index_ptr(neighbor_quads, Ptr{p4est_quadrant_t}, 0),
                    )[],
                )
                sc_array_destroy(neighbor_quads)
                sc_array_destroy(neighbor_encs)
                sc_array_destroy(neighbor_qid)
                return Ptr{PS_Data_2D}(pq.p.user_data[])
            end
        else
            return nothing
        end
    end
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

function DVM_4est_new(
    comm::MPI.Comm,
    conn::Ptr{p4est_connectivity},
    min_quads::Integer,
    min_level::Integer,
    uniform::Integer,
    ::Type{T},
    init_fn::Function,
    user_pointer::Ptr,
) where {T}
    GC.@preserve comm conn min_quads min_level uniform user_pointer p4est_new_ext(
        comm,
        conn,
        min_quads,
        min_level,
        uniform,
        sizeof(T),
        @cfunction($init_fn, Cvoid, (Ptr{p4est_t}, p4est_topidx_t, Ptr{p4est_quadrant_t})),
        user_pointer,
    )
end
function DVM_4est_new(comm::MPI.Comm, conn::Ptr{p4est_connectivity}, ::Type{T}) where {T}
    GC.@preserve comm conn p4est_new_ext(comm, conn, 1, 0, 1, sizeof(T), C_NULL, C_NULL)
end
function DVM_4est_new(
    comm::MPI.Comm,
    conn::Ptr{p4est_connectivity},
    ::Type{T},
    init_fn::Function,
    user_pointer::Ptr,
) where {T}
    GC.@preserve comm conn user_pointer p4est_new_ext(
        comm,
        conn,
        1,
        0,
        1,
        sizeof(T),
        @cfunction($init_fn, Cvoid, (Ptr{p4est_t}, p4est_topidx_t, Ptr{p4est_quadrant_t})),
        user_pointer,
    )
end
function DVM_4est_new(
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
function DVM_4est_new(
    comm::MPI.Comm,
    conn::Ptr{p4est_connectivity},
    ::Type{T},
    init_fn::Function,
) where {T}
    GC.@preserve comm conn user_pointer init_fn p4est_new_ext(
        comm,
        conn,
        1,
        0,
        1,
        sizeof(T),
        @cfunction($init_fn, Cvoid, (Ptr{p4est_t}, p4est_topidx_t, Ptr{p4est_quadrant_t})),
        C_NULL,
    )
end
function DVM_quad_init(
    forest::Ptr{p4est_t},
    which_tree::p4est_topidx_t,
    quadrant::Ptr{p4est_quadrant_t},
    ::Type{T},
    kernel::Function,
) where {T}
    GC.@preserve forest which_tree quadrant begin
        fp = PointerWrapper(forest)
        qp = PointerWrapper(quadrant)
        dp = PointerWrapper(T, qp.p.user_data[])
        GC.@preserve fp qp dp kernel(fp, which_tree, qp, dp)
    end
    return nothing
end
function DVM_4est_volume_iterate(
    forest::Ptr{p4est_t},
    ghost::Ptr{p4est_ghost_t},
    user_data::Ptr{Nothing},
    volume_iter_fn::Function,
)
    GC.@preserve forest ghost user_data volume_iter_fn p4est_iterate(
        forest,
        ghost,
        user_data,
        @cfunction($volume_iter_fn, Cvoid, (Ptr{p4est_iter_volume_info}, Ptr{Nothing})),
        C_NULL,
        C_NULL,
    )
end
function DVM_4est_volume_iterate(
    forest::Ptr{p4est_t},
    user_data::Ptr{Nothing},
    volume_iter_fn::Function,
)
    GC.@preserve forest user_data volume_iter_fn p4est_iterate(
        forest,
        C_NULL,
        user_data,
        @cfunction($volume_iter_fn, Cvoid, (Ptr{p4est_iter_volume_info}, Ptr{Nothing})),
        C_NULL,
        C_NULL,
    )
end
function DVM_volume_iterate(
    info::Ptr{p4est_iter_volume_info},
    data::Ptr{Nothing},
    ::Type{T},
    kernel::Function,
) where {T}
    GC.@preserve info data begin
        ip = PointerWrapper(info)
        dp = PointerWrapper(T, ip.quad.p.user_data[])
        GC.@preserve ip dp kernel(ip, data, dp)
    end
    return nothing
end
function DVM_4est_face_iterate(
    forest::Ptr{p4est_t},
    ghost::Ptr{p4est_ghost_t},
    user_data::Ptr{Nothing},
    face_iter_fn::Ptr{Nothing},
)
    GC.@preserve forest ghost user_data face_iter_fn p4est_iterate(
        forest,
        ghost,
        user_data,
        C_NULL,
        face_iter_fn,
        C_NULL,
    )
end
function DVM_4est_face_iterate(
    forest::Ptr{p4est_t},
    user_data::Ptr{Nothing},
    face_iter_fn::Ptr{Nothing},
)
    GC.@preserve forest user_data face_iter_fn p4est_iterate(
        forest,
        C_NULL,
        user_data,
        C_NULL,
        face_iter_fn,
        C_NULL,
    )
end
function DVM_face_iterate(
    info::Ptr{p4est_iter_face_info},
    data::Ptr{Nothing},
    kernel::Function,
)
    GC.@preserve info data kernel begin
        ip = PointerWrapper(info)
        GC.@preserve ip kernel(ip, data)
    end
    return nothing
end
function global_quadid(ip::PointerWrapper{p4est_iter_volume_info_t})
    tp = get_index_pw(ip.p4est.trees, p4est_tree_t, ip.treeid[])
    global_id = tp.quadrants_offset[] + ip.quadid[]
end

function Base.unsafe_load(pw::PointerWrapper{NTuple{N,T}}, i::Integer = 1) where {N,T}
    Base.unsafe_wrap(Vector{T}, Ptr{T}(pointer(pw)), N)
end

function DVM_refine(
    p4est::Ptr{p4est_t},
    flag_fn::Function,
    init_fn::Function,
    replace_fn::Function,
    recursive::Integer = 0,
    maxlevel::Integer = 1,
)
    GC.@preserve p4est flag_fn init_fn replace_fn recursive maxlevel p4est_refine_ext(
        p4est,
        recursive,
        maxlevel,
        @cfunction($flag_fn, Cint, (Ptr{p4est_t}, p4est_topidx_t, Ptr{p4est_quadrant_t})),
        @cfunction($init_fn, Cvoid, (Ptr{p4est_t}, p4est_topidx_t, Ptr{p4est_quadrant_t})),
        @cfunction(
            $replace_fn,
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
    )
end

function DVM_refine_flag(
    forest::Ptr{p4est_t},
    which_tree,
    quadrant,
    ::Type{T},
    kernel::Function,
) where {T}
    GC.@preserve forest which_tree quadrant begin
        qp = PointerWrapper(quadrant)
        dp = PointerWrapper(T, qp.p.user_data[])
        return Cint(kernel(forest, which_tree, qp, dp))
    end
end
function DVM_4est_new(
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
function DVM_4est_volume_iterate(
    forest::Ptr{p8est_t},
    user_data::Ptr{Nothing},
    volume_iter_fn::Function,
)
    GC.@preserve forest user_data volume_iter_fn p8est_iterate(
        forest,
        C_NULL,
        user_data,
        @cfunction($volume_iter_fn, Cvoid, (Ptr{p8est_iter_volume_info}, Ptr{Nothing})),
        C_NULL,
        C_NULL,
        C_NULL
    )
end
function DVM_volume_iterate(
    info::Ptr{p8est_iter_volume_info},
    data::Ptr{Nothing},
    ::Type{T},
    kernel::Function,
) where {T}
    GC.@preserve info data begin
        ip = PointerWrapper(info)
        dp = PointerWrapper(T, ip.quad.p.user_data[])
        GC.@preserve ip dp kernel(ip, data, dp)
    end
    return nothing
end
function get_midpoint(forest::Ptr{p8est_t}, treeid, quadrant::Ptr{p8est_quadrant_t})
    fp = PointerWrapper(forest)
    qp = PointerWrapper(quadrant)
    get_midpoint(fp, treeid, qp)
end
function get_midpoint(
    fp::PointerWrapper{p8est_t},
    treeid,
    qp::PointerWrapper{p8est_quadrant_t},
)
    midpoint = zeros(3)
    GC.@preserve midpoint fp qp treeid begin
        halflength = P8EST_QUADRANT_LEN(qp.level[]) / 2
        x = qp.x[] + halflength
        y = qp.y[] + halflength
        z = qp.z[] + halflength
        p8est_qcoord_to_vertex(pointer(fp.connectivity), treeid, x, y, z, midpoint)
    end
    midpoint
end