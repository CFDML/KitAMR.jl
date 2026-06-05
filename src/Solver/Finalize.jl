"""
$(TYPEDSIGNATURES)
Update residuals.
"""
function residual_check!(ps_data::PS_Data, prim::Vector{Float64}, global_data::Global_Data)
    Res = global_data.status.residual
    Res.step%RES_CHECK_INTERVAL!=0&&(return nothing)
    @. Res.sumRes+=(prim-ps_data.prim) .^ 2
    @. Res.sumAvg+=abs(prim)
    return nothing
end
function residual_comm!(global_data::Global_Data)
    Res = global_data.status.residual
    fp = PointerWrapper(global_data.forest.p4est)
    N = fp.global_num_quadrants[]
    Res.step%RES_CHECK_INTERVAL!=0&&(return nothing)
    MPI.Reduce!(Res.sumRes, (x, y)->x .+ y, 0, MPI.COMM_WORLD)
    MPI.Reduce!(Res.sumAvg, (x, y)->x .+ y, 0, MPI.COMM_WORLD)
    if MPI.Comm_rank(MPI.COMM_WORLD)==0
        @. Res.residual=sqrt(Res.sumRes*N)/(Res.sumAvg+EPS)
    end
    MPI.Bcast!(Res.residual, 0, MPI.COMM_WORLD)
    Res.sumRes.=0.0;
    Res.sumAvg.=0.0
end
"""
$(SIGNATURES)
Check whether the `residual` and `redundant_step` in [`Status`](@ref) satisfies the convergence criterion.
"""
function check_for_convergence(amr::KitAMR_Data)
    maximum(amr.global_data.status.residual.residual)<TOLERANCE&&(
        amr.global_data.status.residual.redundant_step+=1
    )
    return amr.global_data.status.residual.redundant_step>REDUNDANT_STEPS_NUM
end


function finalize_ghost!(ghost_exchange::Ghost_Exchange)
    id = P4est.package_id()
    sc_free(id, ghost_exchange.ghost_datas)
    sc_free(id, ghost_exchange.ghost_slopes)
    sc_free(id, ghost_exchange.ghost_structures)
    for i in eachindex(ghost_exchange.mirror_data_pointers)
        sc_free(id, ghost_exchange.mirror_data_pointers[i])
        sc_free(id, ghost_exchange.mirror_slope_pointers[i])
        sc_free(id, ghost_exchange.mirror_structure_pointers[i])
    end
end
function finalize_p4est!(ps4est::Ptr{p4est_t}, amr::KitAMR_Data)
    global_data = amr.global_data
    p4est_ghost_destroy(global_data.forest.ghost)
    p4est_mesh_destroy(global_data.forest.mesh)
    pp = PointerWrapper(ps4est)
    p4est_connectivity_destroy(pointer(pp.connectivity))
    p4est_destroy(ps4est)
end

"""
$(TYPEDSIGNATURES)
"""
function finalize!(ps4est::Ptr{p4est_t}, amr::KitAMR_Data)
    finalize_ghost!(amr.ghost.ghost_exchange)
    finalize_p4est!(ps4est, amr)
    return nothing
end

function finalize_p4est!(ps4est::Ptr{p8est_t}, amr::KitAMR_Data)
    global_data = amr.global_data
    p8est_ghost_destroy(global_data.forest.ghost)
    p8est_mesh_destroy(global_data.forest.mesh)
    pp = PointerWrapper(ps4est)
    p8est_connectivity_destroy(pointer(pp.connectivity))
    p8est_destroy(ps4est)
end

"""
$(TYPEDSIGNATURES)
"""
function finalize!(ps4est::Ptr{p8est_t}, amr::KitAMR_Data)
    finalize_ghost!(amr.ghost.ghost_exchange)
    finalize_p4est!(ps4est, amr)
    return nothing
end
