# access_neighbor_w(data::PointerWrapper) = data.w[]
# access_neighbor_w(data::PsData_2D) = data.w
function access_neighbor(
    p4est::Ptr{p4est_t},
    quadid::p4est_locidx_t,
    kinfo::KInfo{DIM,NDF},
    ghost_wrap::Array{AbstractGhostPsData},
    dir::Integer,
) where{DIM,NDF}
    neighbor_quads = sc_array_new(sizeof(Ptr{p4est_quadrant_t}))
    neighbor_encs = sc_array_new(sizeof(Cint))
    neighbor_qid = sc_array_new(sizeof(Cint))
    ghost = kinfo.forest.ghost
    mesh = kinfo.forest.mesh
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
        neighbor = similar(encs, NeighborQuad{DIM,NDF})
        state = 0
        (length(encs) == 0 || encs === nothing) && return [nothing], state # boundary:: type = nothing, state = 0
        for i in eachindex(encs)
            if encs[i] > 0 && encs[i] < 9
                pq = PointerWrapper(
                    iPointerWrapper(neighbor_quads, Ptr{p4est_quadrant_t}, i - 1)[],
                )
                dp = PointerWrapper(P4estPsData, pq.p.user_data[])
                neighbor[i] = unsafe_pointer_to_objref(pointer(dp.ps_data))
                state = 1
            elseif encs[i] > -9 && encs[i] < 0
                neighbor[i] = ghost_wrap[qid[i]+1]
                state = 1
            elseif encs[i] > 8 && encs[i] < 25
                pq = PointerWrapper(
                    iPointerWrapper(neighbor_quads, Ptr{p4est_quadrant_t}, i - 1)[],
                )
                dp = PointerWrapper(P4estPsData, pq.p.user_data[])
                neighbor[i] = unsafe_pointer_to_objref(pointer(dp.ps_data))
                state = -1
            elseif encs[i] > 24
                pq = PointerWrapper(
                    iPointerWrapper(neighbor_quads, Ptr{p4est_quadrant_t}, i - 1)[],
                )
                dp = PointerWrapper(P4estPsData, pq.p.user_data[])
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
    kinfo::KInfo{DIM,NDF},
    ghost_wrap::Array{AbstractGhostPsData},
    dir::Integer,
) where{DIM,NDF}
    neighbor_quads = sc_array_new(sizeof(Ptr{p8est_quadrant_t}))
    neighbor_encs = sc_array_new(sizeof(Cint))
    neighbor_qid = sc_array_new(sizeof(Cint))
    ghost = kinfo.forest.ghost
    mesh = kinfo.forest.mesh
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
        neighbor = similar(encs, NeighborQuad{DIM,NDF})
        state = 0
        (length(encs) == 0 || encs === nothing) && return [nothing], state # boundary:: type = nothing, state = 0
        for i in eachindex(encs)
            if encs[i] > 0 && encs[i] < 25
                pq = PointerWrapper(
                    iPointerWrapper(neighbor_quads, Ptr{p8est_quadrant_t}, i - 1)[],
                )
                dp = PointerWrapper(P4estPsData, pq.p.user_data[])
                neighbor[i] = unsafe_pointer_to_objref(pointer(dp.ps_data))
                state = 1
            elseif encs[i] > -25 && encs[i] < 0
                neighbor[i] = ghost_wrap[qid[i]+1]
                state = 1
            elseif encs[i] > 24 && encs[i] < 121
                pq = PointerWrapper(
                    iPointerWrapper(neighbor_quads, Ptr{p8est_quadrant_t}, i - 1)[],
                )
                dp = PointerWrapper(P4estPsData, pq.p.user_data[])
                neighbor[i] = unsafe_pointer_to_objref(pointer(dp.ps_data))
                state = -1
            elseif encs[i] > 120
                pq = PointerWrapper(
                    iPointerWrapper(neighbor_quads, Ptr{p8est_quadrant_t}, i - 1)[],
                )
                dp = PointerWrapper(P4estPsData, pq.p.user_data[])
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

"""
$(TYPEDSIGNATURES)
"""
function initialize_neighbor_data!(p4est::P_pxest_t, ka::KA)
    p_ka = pointer_from_objref(ka)
    GC.@preserve  ka AMR_volume_iterate(p4est;ghost = ka.kinfo.forest.ghost,user_data = p_ka) do ip,data,dp
        ka = unsafe_pointer_to_objref(data)
        ps_data = unsafe_pointer_to_objref(pointer(dp.ps_data))
        isa(ps_data,InsideSolidData) && return nothing
        face_num = ip isa PointerWrapper{p4est_iter_volume_info_t} ? face_num_2d : face_num_3d
        for i = 1:face_num
            ps_data.neighbor.data[i], ps_data.neighbor.state[i] = access_neighbor(
                pointer(ip.p4est),
                local_quadid(ip),
                ka.kinfo,
                ka.kdata.ghost.ghost_wrap,
                i - 1,
            )
        end
        corner_range = ip isa PointerWrapper{p4est_iter_volume_info_t} ? (4:7) : (6:25)
        if ps_data.bound_enc<0
            for i in corner_range
                nb,state = access_neighbor(
                    pointer(ip.p4est),
                    local_quadid(ip),
                    ka.kinfo,
                    ka.kdata.ghost.ghost_wrap,
                    i,
                )
                push!(ps_data.neighbor.data,nb);push!(ps_data.neighbor.state,state)
            end
        end
        return nothing
    end
end

#=
face_micro_nums: 1->4时记1，只需在state=-1时+0.25即可
=#
function update_neighbor_kernel!(p4est::P_pxest_t, ka::KA)
    p_ka = pointer_from_objref(ka)
    MPI.Barrier(MPI.COMM_WORLD)
    GC.@preserve  ka AMR_volume_iterate(p4est;ghost = ka.kinfo.forest.ghost,user_data = p_ka) do ip,data,dp
        ka = unsafe_pointer_to_objref(data)
        ps_data = unsafe_pointer_to_objref(pointer(dp.ps_data))
        isa(ps_data,InsideSolidData) && return nothing
        is2d = ip isa PointerWrapper{p4est_iter_volume_info_t}
        face_num = is2d ? face_num_2d : face_num_3d
        extra_range = is2d ? (4:7) : (6:25)
        extra_thresh = is2d ? 4 : 6
        for i = 1:face_num
            ps_data.neighbor.data[i], ps_data.neighbor.state[i] = access_neighbor(
                pointer(ip.p4est),
                local_quadid(ip),
                ka.kinfo,
                ka.kdata.ghost.ghost_wrap,
                i - 1,
            )
        end
        if length(ps_data.neighbor.data) > extra_thresh
            for i in extra_range
                ps_data.neighbor.data[i+1],ps_data.neighbor.state[i+1] = access_neighbor(
                    pointer(ip.p4est),
                    local_quadid(ip),
                    ka.kinfo,
                    ka.kdata.ghost.ghost_wrap,
                    i,
                )
            end
        elseif ps_data.bound_enc<0
            for i in extra_range
                nb,state = access_neighbor(
                    pointer(ip.p4est),
                    local_quadid(ip),
                    ka.kinfo,
                    ka.kdata.ghost.ghost_wrap,
                    i,
                )
                push!(ps_data.neighbor.data,nb);push!(ps_data.neighbor.state,state)
            end
        end
        return nothing
    end
end
function update_neighbor!(p4est::Ptr{p4est_t}, ka::KA)
    kinfo = ka.kinfo
    p4est_mesh_destroy(kinfo.forest.mesh)
    kinfo.forest.mesh =
        p4est_mesh_new_ext(p4est, kinfo.forest.ghost, 1, 1,  P4EST_CONNECT_FULL)
    update_neighbor_kernel!(p4est, ka)
end
function update_neighbor!(p4est::Ptr{p8est_t}, ka::KA)
    kinfo = ka.kinfo
    p8est_mesh_destroy(kinfo.forest.mesh)
    kinfo.forest.mesh =
        p8est_mesh_new_ext(p4est, kinfo.forest.ghost, 1, 1, P8EST_CONNECT_FULL)
    update_neighbor_kernel!(p4est, ka)
end
