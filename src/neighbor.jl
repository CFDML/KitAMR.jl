# access_neighbor_w(data::PointerWrapper) = data.w[]
# access_neighbor_w(data::PS_Data_2D) = data.w
function access_neighbor(
    p4est::Ptr{p4est_t},
    quadid::p4est_locidx_t,
    global_data::Global_Data{DIM,NDF},
    ghost_wrap::Array{Ghost_PS_Data},
    dir::Integer,
) where{DIM,NDF}
    neighbor_quads = sc_array_new(sizeof(Ptr{p4est_quadrant_t}))
    neighbor_encs = sc_array_new(sizeof(Cint))
    neighbor_qid = sc_array_new(sizeof(Cint))
    ghost = global_data.forest.ghost
    mesh = global_data.forest.mesh
    GC.@preserve p4est ghost mesh quadid dir neighbor_encs neighbor_qid neighbor_quads begin
        p4est_mesh_get_neighbors(
            p4est,
            ghost,
            mesh,
            quadid,
            dir,
            neighbor_quads,
            neighbor_encs,
            neighbor_qid,
        )
        qid = unsafe_wrap_sc(Cint, neighbor_qid)
        encs = unsafe_wrap_sc(Cint, neighbor_encs)
        neighbor = similar(encs, Neighbor_Quad{DIM,NDF})
        state = 0
        (length(encs) == 0 || encs === nothing) && return [nothing], state # boundary:: type = nothing, state = 0
        for i in eachindex(encs)
            if encs[i] > 0 && encs[i] < 9
                pq = PointerWrapper(
                    iPointerWrapper(neighbor_quads, Ptr{p4est_quadrant_t}, i - 1)[],
                )
                dp = PointerWrapper(P4est_PS_Data, pq.p.user_data[])
                neighbor[i] = unsafe_pointer_to_objref(pointer(dp.ps_data))
                state = 1
            elseif encs[i] > -9 && encs[i] < 0
                neighbor[i] = ghost_wrap[qid[i]+1]
                state = 1
            elseif encs[i] > 8 && encs[i] < 25
                pq = PointerWrapper(
                    iPointerWrapper(neighbor_quads, Ptr{p4est_quadrant_t}, i - 1)[],
                )
                dp = PointerWrapper(P4est_PS_Data, pq.p.user_data[])
                neighbor[i] = unsafe_pointer_to_objref(pointer(dp.ps_data))
                state = -1
            elseif encs[i] > 24
                pq = PointerWrapper(
                    iPointerWrapper(neighbor_quads, Ptr{p4est_quadrant_t}, i - 1)[],
                )
                dp = PointerWrapper(P4est_PS_Data, pq.p.user_data[])
                neighbor[i] = unsafe_pointer_to_objref(pointer(dp.ps_data))
                state = 2^(DIM - 1)
            elseif encs[i] > -25 && encs[i] < -8
                neighbor[i] = ghost_wrap[qid[i]+1]
                state = -1
            elseif encs[i] < -24
                neighbor[i] = ghost_wrap[qid[i]+1]
                state = 2^(DIM - 1)
            end
        end
    end
    sc_array_destroy(neighbor_quads)
    sc_array_destroy(neighbor_encs)
    sc_array_destroy(neighbor_qid)
    return (neighbor, state)
end
function access_neighbor(
    p4est::Ptr{p8est_t},
    quadid::p4est_locidx_t,
    global_data::Global_Data{DIM,NDF},
    ghost_wrap::Array{Ghost_PS_Data},
    dir::Integer,
) where{DIM,NDF}
    neighbor_quads = sc_array_new(sizeof(Ptr{p8est_quadrant_t}))
    neighbor_encs = sc_array_new(sizeof(Cint))
    neighbor_qid = sc_array_new(sizeof(Cint))
    ghost = global_data.forest.ghost
    mesh = global_data.forest.mesh
    GC.@preserve p4est ghost mesh quadid dir neighbor_encs neighbor_qid neighbor_quads begin
        p8est_mesh_get_neighbors(
            p4est,
            ghost,
            mesh,
            quadid,
            dir,
            neighbor_quads,
            neighbor_encs,
            neighbor_qid,
        )
        qid = unsafe_wrap_sc(Cint, neighbor_qid)
        encs = unsafe_wrap_sc(Cint, neighbor_encs)
        neighbor = similar(encs, Neighbor_Quad{DIM,NDF})
        state = 0
        (length(encs) == 0 || encs === nothing) && return [nothing], state # boundary:: type = nothing, state = 0
        for i in eachindex(encs)
            if encs[i] > 0 && encs[i] < 25
                pq = PointerWrapper(
                    iPointerWrapper(neighbor_quads, Ptr{p8est_quadrant_t}, i - 1)[],
                )
                dp = PointerWrapper(P4est_PS_Data, pq.p.user_data[])
                neighbor[i] = unsafe_pointer_to_objref(pointer(dp.ps_data))
                state = 1
            elseif encs[i] > -25 && encs[i] < 0
                neighbor[i] = ghost_wrap[qid[i]+1]
                state = 1
            elseif encs[i] > 24 && encs[i] < 121
                pq = PointerWrapper(
                    iPointerWrapper(neighbor_quads, Ptr{p8est_quadrant_t}, i - 1)[],
                )
                dp = PointerWrapper(P4est_PS_Data, pq.p.user_data[])
                neighbor[i] = unsafe_pointer_to_objref(pointer(dp.ps_data))
                state = -1
            elseif encs[i] > 120
                pq = PointerWrapper(
                    iPointerWrapper(neighbor_quads, Ptr{p8est_quadrant_t}, i - 1)[],
                )
                dp = PointerWrapper(P4est_PS_Data, pq.p.user_data[])
                neighbor[i] = unsafe_pointer_to_objref(pointer(dp.ps_data))
                state = 2^(DIM - 1)
            elseif encs[i] > -121 && encs[i] < -24
                neighbor[i] = ghost_wrap[qid[i]+1]
                state = -1
            elseif encs[i] < -120
                neighbor[i] = ghost_wrap[qid[i]+1]
                state = 2^(DIM - 1)
            end
        end
    end
    sc_array_destroy(neighbor_quads)
    sc_array_destroy(neighbor_encs)
    sc_array_destroy(neighbor_qid)
    return (neighbor, state)
end
function initialize_neighbor_data!(ip::PointerWrapper{p4est_iter_volume_info_t}, data, dp)
    amr = unsafe_pointer_to_objref(data)
    ps_data = unsafe_pointer_to_objref(pointer(dp.ps_data))
    for i = 1:face_num_2d
        ps_data.neighbor.data[i], ps_data.neighbor.state[i] = access_neighbor(
            pointer(ip.p4est),
            global_quadid(ip),
            amr.global_data,
            amr.ghost.ghost_wrap,
            i - 1,
        )
    end
end
function initialize_neighbor_data!(ip::PointerWrapper{p8est_iter_volume_info_t}, data, dp)
    amr = unsafe_pointer_to_objref(data)
    ps_data = unsafe_pointer_to_objref(pointer(dp.ps_data))
    for i = 1:face_num_3d
        ps_data.neighbor.data[i], ps_data.neighbor.state[i] = access_neighbor(
            pointer(ip.p4est),
            global_quadid(ip),
            amr.global_data,
            amr.ghost.ghost_wrap,
            i - 1,
        )
    end
end
function initialize_neighbor_data!(info, data)
    DVM_volume_iterate(info, data, P4est_PS_Data, initialize_neighbor_data!)
end
function initialize_neighbor_data!(ps4est::P_pxest_t, amr::AMR)
    p_amr = pointer_from_objref(amr)
    GC.@preserve amr DVM_4est_volume_iterate(
        ps4est,
        amr.global_data.forest.ghost,
        p_amr,
        initialize_neighbor_data!,
    )
end
#=
face_micro_nums: 1->4时记1，只需在state=-1时+0.25即可
=#
function update_neighbor_kernel!(ip::PointerWrapper{p4est_iter_volume_info_t}, data, dp)
    amr = unsafe_pointer_to_objref(data)
    ps_data = unsafe_pointer_to_objref(pointer(dp.ps_data))
    for i = 1:face_num_2d
        ps_data.neighbor.data[i], ps_data.neighbor.state[i] = access_neighbor(
            pointer(ip.p4est),
            global_quadid(ip),
            amr.global_data,
            amr.ghost.ghost_wrap,
            i - 1,
        )
    end
end
function update_neighbor_kernel!(ip::PointerWrapper{p8est_iter_volume_info_t}, data, dp)
    amr = unsafe_pointer_to_objref(data)
    ps_data = unsafe_pointer_to_objref(pointer(dp.ps_data))
    for i = 1:face_num_3d
        ps_data.neighbor.data[i], ps_data.neighbor.state[i] = access_neighbor(
            pointer(ip.p4est),
            global_quadid(ip),
            amr.global_data,
            amr.ghost.ghost_wrap,
            i - 1,
        )
    end
end
function update_neighbor_kernel!(info, data)
    DVM_volume_iterate(info, data, P4est_PS_Data, update_neighbor_kernel!)
end
function update_neighbor_kernel!(ps4est::P_pxest_t, amr::AMR)
    p_amr = pointer_from_objref(amr)
    GC.@preserve amr DVM_4est_volume_iterate(
        ps4est,
        amr.global_data.forest.ghost,
        p_amr,
        update_neighbor_kernel!,
    )
end
function update_neighbor!(p4est::Ptr{p4est_t}, amr::AMR)
    global_data = amr.global_data
    p4est_mesh_destroy(global_data.forest.mesh)
    global_data.forest.mesh =
        p4est_mesh_new_ext(p4est, global_data.forest.ghost, 1, 1, P4EST_CONNECT_FACE)
    update_neighbor_kernel!(p4est, amr)
end
function update_neighbor!(p4est::Ptr{p8est_t}, amr::AMR)
    global_data = amr.global_data
    p8est_mesh_destroy(global_data.forest.mesh)
    global_data.forest.mesh =
        p8est_mesh_new_ext(p4est, global_data.forest.ghost, 1, 1, P8EST_CONNECT_FACE)
    update_neighbor_kernel!(p4est, amr)
end
