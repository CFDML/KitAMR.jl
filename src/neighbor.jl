# access_neighbor_w(data::PointerWrapper) = data.w[]
# access_neighbor_w(data::PS_Data) = data.w
function access_neighbor(
    p4est::Ptr{p4est_t},
    quadid::p4est_locidx_t,
    global_data::Global_Data,
    ghost_wrap::Array{Ghost_PS_Data},
    dir::Integer,
)
    neighbor_quads = sc_array_new(sizeof(Ptr{p4est_quadrant_t}))
    neighbor_encs = sc_array_new(sizeof(Cint))
    neighbor_qid = sc_array_new(sizeof(Cint))
    ghost = global_data.ghost
    mesh = global_data.mesh
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
        neighbor = similar(encs, Union{AbstractPsData,Nothing})
        state = 0
        (length(encs) == 0 || encs === nothing) && return [nothing], state # boundary:: type = nothing, state = 1
        for i in eachindex(encs)
            if encs[i] > 0 && encs[i] < 9
                pq = PointerWrapper(
                    get_index_pw(neighbor_quads, Ptr{p4est_quadrant_t}, i - 1)[],
                )
                dp = PointerWrapper(P4est_PS_Data, pq.p.user_data[])
                neighbor[i] = unsafe_pointer_to_objref(pointer(dp.ps_data))
                state = 1
            elseif encs[i] > -9 && encs[i] < 0
                neighbor[i] = ghost_wrap[qid[i]+1]
                state = 1
            elseif encs[i] > 8 && encs[i] < 25
                pq = PointerWrapper(
                    get_index_pw(neighbor_quads, Ptr{p4est_quadrant_t}, i - 1)[],
                )
                dp = PointerWrapper(P4est_PS_Data, pq.p.user_data[])
                neighbor[i] = unsafe_pointer_to_objref(pointer(dp.ps_data))
                state = -1
            elseif encs[i] > 24
                pq = PointerWrapper(
                    get_index_pw(neighbor_quads, Ptr{p4est_quadrant_t}, i - 1)[],
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
function initialize_neighbor_data!(ip::PointerWrapper{p4est_iter_volume_info_t}, data, dp)
    DVM_data = unsafe_pointer_to_objref(data)
    ps_data = unsafe_pointer_to_objref(pointer(dp.ps_data))
    for i = 1:2*DIM
        ps_data.neighbor.data[i], ps_data.neighbor.state[i] = access_neighbor(
            pointer(ip.p4est),
            global_quadid(ip),
            DVM_data.global_data,
            DVM_data.ghost_wrap,
            i - 1,
        )
    end
end
function initialize_neighbor_data!(info, data)
    DVM_volume_iterate(info, data, P4est_PS_Data, initialize_neighbor_data!)
end
function initialize_neighbor_data!(ps4est::Ptr{p4est_t}, DVM_data::DVM_Data)
    p_DVM_data = pointer_from_objref(DVM_data)
    GC.@preserve DVM_data DVM_4est_volume_iterate(
        ps4est,
        DVM_data.global_data.ghost,
        p_DVM_data,
        initialize_neighbor_data!,
    )
end
#=
face_micro_nums: 1->4时记1，只需在state=-1时+0.25即可
=#
function update_neighbor_kernel!(ip::PointerWrapper{p4est_iter_volume_info_t}, data, dp)
    DVM_data = unsafe_pointer_to_objref(data)
    ps_data = unsafe_pointer_to_objref(pointer(dp.ps_data))
    for i = 1:2*DIM
        ps_data.neighbor.data[i], ps_data.neighbor.state[i] = access_neighbor(
            pointer(ip.p4est),
            global_quadid(ip),
            DVM_data.global_data,
            DVM_data.ghost_wrap,
            i - 1,
        )
    end
end
function update_neighbor_kernel!(info, data)
    DVM_volume_iterate(info, data, P4est_PS_Data, update_neighbor_kernel!)
end
function update_neighbor_kernel!(ps4est::Ptr{p4est_t}, DVM_data::DVM_Data)
    p_DVM_data = pointer_from_objref(DVM_data)
    GC.@preserve DVM_data DVM_4est_volume_iterate(
        ps4est,
        DVM_data.global_data.ghost,
        p_DVM_data,
        update_neighbor_kernel!,
    )
end
function update_neighbor!(p4est::Ptr{p4est_t}, DVM_data::DVM_Data)
    global_data = DVM_data.global_data
    p4est_mesh_destroy(global_data.mesh)
    global_data.mesh =
        p4est_mesh_new_ext(p4est, global_data.ghost, 1, 1, P4EST_CONNECT_FACE)
    update_neighbor_kernel!(p4est, DVM_data)
end
