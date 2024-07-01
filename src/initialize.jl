function initialize_faces!(
    ::Val{0},
    ::Val{1},
    ip::PointerWrapper{p4est_iter_face_info_t},
    side,
    DVM_data,
)
    faces = DVM_data.faces
    base_quad =
        pointer(PointerWrapper(P4est_PS_Data, side.is.full.quad.p.user_data[]).ps_data)
    faceid = side.face[] + 1
    push!(faces, Face(unsafe_pointer_to_objref(base_quad), faceid, 0, nothing))
end
function initialize_faces!(
    ::Val{0},
    ::Val{0},
    ip::PointerWrapper{p4est_iter_face_info_t},
    side,
    DVM_data,
)
    faces = DVM_data.faces
    base_quad =
        pointer(PointerWrapper(P4est_PS_Data, side.is.full.quad.p.user_data[]).ps_data)
    faceid = side.face[] + 1
    push!(faces, Face(unsafe_pointer_to_objref(base_quad), faceid, 0, nothing))
end
function initialize_faces!(
    ::Val{1},
    ::Val{1},
    ip::PointerWrapper{p4est_iter_face_info_t},
    side,
    DVM_data,
)
    faces = DVM_data.faces
    side = get_index_pw(ip.sides, p4est_iter_face_side_t, 0)
    is_ghost = Base.unsafe_wrap(
        Vector{Int8},
        Ptr{Int8}(pointer(side.is.hanging.is_ghost)),
        2^(DIM - 1),
    )
    quadid = Base.unsafe_wrap(
        Vector{Int32},
        Ptr{Int32}(pointer(side.is.hanging.quadid)),
        2^(DIM - 1),
    )
    # base_quad = C_NULL;hanging_quad = C_NULL;faceid = 0
    baseid = findfirst(x -> x == 0, is_ghost) - 1
    hangingid = baseid == 0 ? 1 : 0
    qp = PointerWrapper(get_index_pw(side.is.hanging.quad, Ptr{p4est_quadrant_t}, baseid)[])
    base_quad = pointer(PointerWrapper(P4est_PS_Data, qp.p.user_data[]).ps_data)
    faceid = side.face[] + 1
    qp = PointerWrapper(
        get_index_pw(side.is.hanging.quad, Ptr{p4est_quadrant_t}, hangingid)[],
    )
    if is_ghost[hangingid+1] == 0
        hanging_quad = unsafe_pointer_to_objref(
            pointer(PointerWrapper(P4est_PS_Data, qp.p.user_data[]).ps_data),
        )
    else
        hanging_quad = DVM_data.ghost_wrap[quadid[hangingid+1]+1]
    end
    push!(faces, Face(unsafe_pointer_to_objref(base_quad), faceid, 2, hanging_quad))
end
function initialize_faces!(
    ::Val{1},
    ::Val{0},
    ip::PointerWrapper{p4est_iter_face_info_t},
    side,
    DVM_data,
)
    side = get_index_pw(ip.sides, p4est_iter_face_side_t, 1)
    faces = DVM_data.faces
    if side.is_hanging[] == 0
        base_quad =
            pointer(PointerWrapper(P4est_PS_Data, side.is.full.quad.p.user_data[]).ps_data)
        faceid = side.face[] + 1
        push!(faces, Face(unsafe_pointer_to_objref(base_quad), faceid, 0, nothing))
    else
        is_ghost = Base.unsafe_wrap(
            Vector{Int8},
            Ptr{Int8}(pointer(side.is.hanging.is_ghost)),
            2^(DIM - 1),
        )
        quadid = Base.unsafe_wrap(
            Vector{Int32},
            Ptr{Int32}(pointer(side.is.hanging.quadid)),
            2^(DIM - 1),
        )
        # base_quad = C_NULL;hanging_quad = C_NULL;faceid = 0
        baseid = findfirst(x -> x == 0, is_ghost) - 1
        hangingid = baseid == 0 ? 1 : 0
        qp = PointerWrapper(
            get_index_pw(side.is.hanging.quad, Ptr{p4est_quadrant_t}, baseid)[],
        )
        base_quad = pointer(PointerWrapper(P4est_PS_Data, qp.p.user_data[]).ps_data)
        faceid = side.face[] + 1
        qp = PointerWrapper(
            get_index_pw(side.is.hanging.quad, Ptr{p4est_quadrant_t}, hangingid)[],
        )
        if is_ghost[hangingid+1] == 0
            hanging_quad = unsafe_pointer_to_objref(
                pointer(PointerWrapper(P4est_PS_Data, qp.p.user_data[]).ps_data),
            )
        else
            hanging_quad = DVM_data.ghost_wrap[quadid[hangingid+1]+1]
        end
        push!(faces, Face(unsafe_pointer_to_objref(base_quad), faceid, 2, hanging_quad))
    end
end
function initialize_faces!(
    ::Val{2},
    ip::PointerWrapper{p4est_iter_face_info_t},
    data::Ptr{Nothing},
)
    DVM_data = unsafe_pointer_to_objref(data)
    sideid = 0
    side = get_index_pw(ip.sides, p4est_iter_face_side_t, sideid)
    side.is_hanging[] == 1 &&
        (side = get_index_pw(ip.sides, p4est_iter_face_side_t, (sideid += 1)))
    # faceid = side.face[]+1
    # if side.is_hanging[]==0
    #     base_quad = pointer(PointerWrapper(P4est_PS_Data, side.is.full.quad.p.user_data[]).ps_data)
    #     push!(faces,Face(unsafe_pointer_to_objref(base_quad),faceid,0,nothing))
    # else
    #     qp = PointerWrapper(get_index_pw(side.is.hanging.quad, Ptr{p4est_quadrant_t}, 0)[])
    #     base_quad = pointer(PointerWrapper(P4est_PS_Data, qp.p.user_data[]).ps_data)
    #     qp = PointerWrapper(get_index_pw(side.is.hanging.quad,Ptr{p4est_quadrant_t}, 1)[])
    #     hanging_quad = pointer(PointerWrapper(P4est_PS_Data, qp.p.user_data[]).ps_data)
    #     push!(faces,Face(unsafe_pointer_to_objref(base_quad),faceid,2,unsafe_pointer_to_objref(hanging_quad)))
    # end
    initialize_faces!(Val(Int(side.is.full.is_ghost[])), Val(sideid), ip, side, DVM_data)
end
function initialize_faces!(
    ::Val{1},
    ip::PointerWrapper{p4est_iter_face_info_t},
    data::Ptr{Nothing},
)
    DVM_data = unsafe_pointer_to_objref(data)
    faces = DVM_data.faces
    side = get_index_pw(ip.sides, p4est_iter_face_side_t, 0)
    base_quad =
        pointer(PointerWrapper(P4est_PS_Data, side.is.full.quad.p.user_data[]).ps_data)
    faceid = side.face[] + 1
    push!(faces, Face(unsafe_pointer_to_objref(base_quad), faceid, 1, nothing))
end
function initialize_faces!(ip::PointerWrapper{p4est_iter_face_info_t}, data)
    initialize_faces!(Val(Int(ip.sides.elem_count[])), ip, data)
end
function initialize_faces!(info::Ptr{p4est_iter_face_info_t}, data)
    DVM_face_iterate(info, data, initialize_faces!)
end
function initialize_faces!(ps4est, DVM_data)
    global_data = DVM_data.global_data
    p_data = pointer_from_objref(DVM_data)
    c_initialize_faces =
        @cfunction(initialize_faces!, Cvoid, (Ptr{p4est_iter_face_info_t}, Ptr{Nothing}))
    GC.@preserve p_data c_initialize_faces DVM_4est_face_iterate(
        ps4est,
        global_data.ghost,
        p_data,
        c_initialize_faces,
    )
end
function init_p4est_kernal(ip, data, dp)
    global_data, trees = unsafe_pointer_to_objref(data)
    treeid = ip.treeid[] - trees.offset
    ps_data = PS_Data()
    push!(trees.data[treeid], ps_data)
    dp[] = P4est_PS_Data(pointer_from_objref(ps_data))
    # geometry = global_data.geometry
    ic = global_data.ic
    # bc = global_data.bc
    ds, midpoint = get_midpoint_ds(ip.p4est, ip.treeid[], ip.quad)
    ps_data.ds .= ds
    ps_data.midpoint .= midpoint
    ps_data.prim .= ic
    ps_data.w .= get_conserved(ps_data.prim, global_data.gas.γ)
    ps_data.vs_data = init_VS(ps_data.prim, global_data)
end
function init_ps4est(info, data)
    DVM_volume_iterate(info, data, P4est_PS_Data, init_p4est_kernal)
end
function init_PS(comm, global_data)
    init_PS(comm, global_data, global_data.trees_num...)
end
# function test_partition(ip,data,dp)
#     # MPI.Comm_rank(MPI.COMM_WORLD)==1 && @show global_quadid(ip)
#     if global_quadid(ip)==0&&MPI.Comm_rank(MPI.COMM_WORLD)==1
#         @show get_midpoint_ds(ip.p4est, ip.treeid[], ip.quad)
#     end
#     #=
#     get_midpoint_ds(ip.p4est, ip.treeid[], ip.quad) = ([0.0125, 0.012500000000000067], [0.00625, 0.29375000000000007])
#     gfq = [0, 418, 677, 1096]
#     gfq = [0, 365, 730, 1096]
#     =#
# end
# function test_partition(info,data)
#     DVM_volume_iterate(info,data,P4est_PS_Data,test_partition)
# end
function re_init_vs4est!(trees, global_data)
    for i in eachindex(trees.data)
        for j in eachindex(trees.data[i])
            ps_data = trees.data[i][j]
            vs_data = ps_data.vs_data
            vs_data.df .=
                discrete_maxwell(vs_data.midpoint, ps_data.prim, global_data.gas.K)
        end
    end
end
function init_PS(comm::MPI.Comm, global_data::Global_Data, Nx::Integer, Ny::Integer)
    GC.@preserve global_data begin
        connectivity_ps = Cartesian_connectivity(Nx, Ny, global_data.geometry...)
        ps4est = DVM_4est_new(
            comm,
            connectivity_ps.pointer,
            P4est_PS_Data,
            pointer_from_objref(global_data),
        )
        # pp = PointerWrapper(ps4est)
        # if MPI.Comm_rank(comm) == 0
        #     gfq = Base.unsafe_wrap(Vector{Int},pointer(pp.global_first_quadrant),MPI.Comm_size(MPI.COMM_WORLD)+1)
        #     @show gfq
        # end
        pre_PS_refine!(ps4est)
        pre_PS_balance!(ps4est)
        p4est_partition(ps4est, 0, C_NULL)
        # DVM_4est_volume_iterate(ps4est,C_NULL,test_partition)
        # if MPI.Comm_rank(comm) == 0
        #     gfq = Base.unsafe_wrap(Vector{Int},pointer(pp.global_first_quadrant),MPI.Comm_size(MPI.COMM_WORLD)+1)
        #     @show gfq
        # end
        fp = PointerWrapper(ps4est)
        trees_data =
            Vector{Vector{PS_Data}}(undef, fp.last_local_tree[] - fp.first_local_tree[] + 1)
        for i in eachindex(trees_data)
            trees_data[i] = Vector{PS_Data}(undef, 0)
        end
        trees = Trees(trees_data, fp.first_local_tree[] - 1)
        data = [global_data, trees]
        p_data = pointer_from_objref(data)
        GC.@preserve data DVM_4est_volume_iterate(ps4est, p_data, init_ps4est)
        pre_vs_refine!(trees, global_data)
        re_init_vs4est!(trees, global_data)
        return trees, ps4est
    end
end
function init()
    global_data = Global_Data()
    trees, ps4est = init_PS(MPI.COMM_WORLD, global_data, global_data.trees_num...)
    GC.@preserve ps4est ghost_ps = p4est_ghost_new(ps4est, P4EST_CONNECT_FACE)
    GC.@preserve ps4est ghost_ps mesh_ps =
        p4est_mesh_new_ext(ps4est, ghost_ps, 1, 1, P4EST_CONNECT_FACE)
    global_data.ghost = ghost_ps
    global_data.mesh = mesh_ps
    ghost_exchange = initialize_ghost_exchange(ps4est, global_data)
    ghost_wrap = initialize_ghost_wrap(global_data, ghost_exchange)
    data_set = Data_Set()
    DVM_data = DVM_Data(
        global_data,
        ghost_exchange,
        ghost_wrap,
        trees,
        Vector{Face}(undef, 0),
        ps4est,
        data_set,
    )
    initialize_faces!(ps4est, DVM_data)
    PointerWrapper(ps4est).user_pointer = pointer_from_objref(DVM_data)
    initialize_neighbor_data!(ps4est, DVM_data)
    return (ps4est, DVM_data)
end
