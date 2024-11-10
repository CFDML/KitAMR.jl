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
    GC.@preserve comm conn min_quads min_level uniform user_pointer p8est_new_ext(
        comm,
        conn,
        min_quads,
        min_level,
        uniform,
        sizeof(T),
        @cfunction($init_fn, Cvoid, (Ptr{p8est_t}, p4est_topidx_t, Ptr{p8est_quadrant_t})),
        user_pointer,
    )
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
function AMR_4est_new(
    comm::MPI.Comm,
    conn::Ptr{p8est_connectivity},
    ::Type{T},
    init_fn::Function,
    user_pointer::Ptr,
) where {T}
    GC.@preserve comm conn user_pointer p8est_new_ext(
        comm,
        conn,
        1,
        0,
        1,
        sizeof(T),
        @cfunction($init_fn, Cvoid, (Ptr{p8est_t}, p4est_topidx_t, Ptr{p8est_quadrant_t})),
        user_pointer,
    )
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
function AMR_4est_new(
    comm::MPI.Comm,
    conn::Ptr{p8est_connectivity},
    ::Type{T},
    init_fn::Function,
) where {T}
    GC.@preserve comm conn user_pointer init_fn p8est_new_ext(
        comm,
        conn,
        1,
        0,
        1,
        sizeof(T),
        @cfunction($init_fn, Cvoid, (Ptr{p8est_t}, p4est_topidx_t, Ptr{p8est_quadrant_t})),
        C_NULL,
    )
end

function AMR_quad_init(
    forest::P_pxest_t,
    which_tree::p4est_topidx_t,
    quadrant::P_pxest_quadrant_t,
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

function AMR_4est_volume_iterate(
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
function AMR_4est_volume_iterate(
    forest::Ptr{p8est_t},
    ghost::Ptr{p8est_ghost_t},
    user_data::Ptr{Nothing},
    volume_iter_fn::Function,
)
    GC.@preserve forest ghost user_data volume_iter_fn p8est_iterate(
        forest,
        ghost,
        user_data,
        @cfunction($volume_iter_fn, Cvoid, (Ptr{p8est_iter_volume_info}, Ptr{Nothing})),
        C_NULL,
        C_NULL,
        C_NULL
    )
end

function AMR_4est_volume_iterate(
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
function AMR_4est_volume_iterate(
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

function AMR_volume_iterate(
    info::P_pxest_iter_volume_info_t,
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

function AMR_4est_face_iterate(
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
function AMR_4est_face_iterate(
    forest::Ptr{p8est_t},
    ghost::Ptr{p8est_ghost_t},
    user_data::Ptr{Nothing},
    face_iter_fn::Ptr{Nothing},
)
    GC.@preserve forest ghost user_data face_iter_fn p8est_iterate(
        forest,
        ghost,
        user_data,
        C_NULL,
        face_iter_fn,
        C_NULL,
        C_NULL
    )
end

function AMR_4est_face_iterate(
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
function AMR_4est_face_iterate(
    forest::Ptr{p8est_t},
    user_data::Ptr{Nothing},
    face_iter_fn::Ptr{Nothing},
)
    GC.@preserve forest user_data face_iter_fn p8est_iterate(
        forest,
        C_NULL,
        user_data,
        C_NULL,
        face_iter_fn,
        C_NULL,
        C_NULL
    )
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

function AMR_ghost_new(p4est::Ptr{p4est_t})
    GC.@preserve p4est p4est_ghost_new(p4est, P4EST_CONNECT_FACE)
end
function AMR_ghost_new(p4est::Ptr{p8est_t})
    GC.@preserve p4est p8est_ghost_new(p4est, P8EST_CONNECT_FACE)
end

function AMR_partition(p4est::Ptr{p4est_t})
    p4est_partition(p4est, 0,C_NULL)
end
function AMR_partition(p4est::Ptr{p8est_t})
    p8est_partition(p4est, 0,C_NULL)
end

function AMR_mesh_new(p4est::Ptr{p4est_t},ghost::Ptr{p4est_ghost_t})
    GC.@preserve p4est ghost p4est_mesh_new_ext(p4est,ghost,1,1,P4EST_CONNECT_FACE)
end
function AMR_mesh_new(p4est::Ptr{p8est_t},ghost::Ptr{p8est_ghost_t})
    GC.@preserve p4est ghost p8est_mesh_new_ext(p4est,ghost,1,1,P8EST_CONNECT_FACE)
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

function AMR_volume_iterate(f::Function,forest::Ptr{p4est_t};ghost=C_NULL,user_data=C_NULL,data_type = P4est_PS_Data)
    function iter_fn(info,data)
        GC.@preserve info data begin
            ip = PointerWrapper(info)
            dp = PointerWrapper(data_type, ip.quad.p.user_data[])
            GC.@preserve ip f(ip, data, dp)
        end
    end
    GC.@preserve forest ghost user_data iter_fn p4est_iterate(
        forest,
        ghost,
        user_data,
        C_NULL,
        @cfunction($iter_fn,Cvoid, (Ptr{p4est_iter_volume_info}, Ptr{Nothing})),
        C_NULL,
    )
end
function AMR_corner_iterate(f::Function,forest::Ptr{p4est_t};ghost=C_NULL,user_data=C_NULL)
    function iter_fn(info,data)
        GC.@preserve info data begin
            ip = PointerWrapper(info)
            GC.@preserve ip f(ip, data)
        end
    end
    GC.@preserve forest ghost user_data iter_fn p4est_iterate(
        forest,
        ghost,
        user_data,
        C_NULL,
        C_NULL,
        @cfunction($iter_fn,Cvoid, (Ptr{p4est_iter_corner_info}, Ptr{Nothing})),
    )
end
function AMR_face_iterate(f::Function,forest::Ptr{p4est_t};ghost=C_NULL,user_data=C_NULL)
    function iter_fn(info,data)
        GC.@preserve info data begin
            ip = PointerWrapper(info)
            GC.@preserve ip f(ip, data)
        end
    end
    GC.@preserve forest ghost user_data iter_fn p4est_iterate(
        forest,
        ghost,
        user_data,
        C_NULL,
        @cfunction($iter_fn,Cvoid, (Ptr{p4est_iter_face_info}, Ptr{Nothing})),
        C_NULL,
    )
end