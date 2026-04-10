"""
$(TYPEDSIGNATURES)
Update residuals.
"""
function residual_check!(ps_data::PsData,prim::Vector{Float64},kinfo::KInfo)
    Res = kinfo.status.residual
    kinfo.status.step%kinfo.config.solver.ST_CHECK_INTERVAL!=0&&(return nothing)
    @. Res.sumRes+=(prim-ps_data.prim).^2
    @. Res.sumAvg+=abs(prim)
    return nothing
end
function residual_comm!(kinfo::KInfo)
    Res = kinfo.status.residual
    fp = PointerWrapper(kinfo.forest.p4est)
    N = fp.global_num_quadrants[]
    kinfo.status.step%kinfo.config.solver.ST_CHECK_INTERVAL!=0&&(return nothing)
    MPI.Reduce!(Res.sumRes,(x,y)->x.+y,0,MPI.COMM_WORLD)
    MPI.Reduce!(Res.sumAvg,(x,y)->x.+y,0,MPI.COMM_WORLD)
    if MPI.Comm_rank(MPI.COMM_WORLD)==0
        @. Res.residual=sqrt(Res.sumRes*N)/(Res.sumAvg+EPS)
    end
    MPI.Bcast!(Res.residual,0,MPI.COMM_WORLD)
    Res.sumRes.=0.;Res.sumAvg.=0.
end
"""
$(SIGNATURES)
Check whether the `residual` and `redundant_step` in [`Status`](@ref) satisfies the convergence criterion.
"""
function check_for_convergence(ka::KA)
    maximum(ka.kinfo.status.residual.residual)<ka.kinfo.config.solver.TOLERANCE&&(ka.kinfo.status.residual.redundant_step+=1)
    return ka.kinfo.status.residual.redundant_step>ka.kinfo.config.solver.REDUNDANT_STEPS_NUM
end


function finalize_ghost!(::GhostBuffer)
    # Julia GC owns all buffers — nothing to free.
    return nothing
end
function finalize_p4est!(p4est::Ptr{p4est_t}, ka::KA)
    kinfo = ka.kinfo
    p4est_ghost_destroy(kinfo.forest.ghost)
    p4est_mesh_destroy(kinfo.forest.mesh)
    pp = PointerWrapper(p4est)
    p4est_connectivity_destroy(pointer(pp.connectivity))
    p4est_destroy(p4est)
end

"""
$(TYPEDSIGNATURES)
"""
function finalize!(p4est::Ptr{p4est_t}, ka::KA)
    finalize_ghost!(ka.kdata.ghost.ghost_buffer)
    finalize_p4est!(p4est, ka)
    return nothing
end

function finalize_p4est!(p4est::Ptr{p8est_t}, ka::KA)
    kinfo = ka.kinfo
    p8est_ghost_destroy(kinfo.forest.ghost)
    p8est_mesh_destroy(kinfo.forest.mesh)
    pp = PointerWrapper(p4est)
    p8est_connectivity_destroy(pointer(pp.connectivity))
    p8est_destroy(p4est)
end

"""
$(TYPEDSIGNATURES)
"""
function finalize!(p4est::Ptr{p8est_t}, ka::KA)
    finalize_ghost!(ka.kdata.ghost.ghost_buffer)
    finalize_p4est!(p4est, ka)
    return nothing
end
