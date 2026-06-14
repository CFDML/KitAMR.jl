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
    MPI.Reduce!(Res.sumRes,+,0,MPI.COMM_WORLD)
    MPI.Reduce!(Res.sumAvg,+,0,MPI.COMM_WORLD)
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
    conn = pointer(PointerWrapper(p4est).connectivity)
    if kinfo.forest.mesh != Ptr{p4est_mesh_t}(C_NULL)
        p4est_mesh_destroy(kinfo.forest.mesh)
        kinfo.forest.mesh = Ptr{p4est_mesh_t}(C_NULL)
    end
    if kinfo.forest.ghost != Ptr{p4est_ghost_t}(C_NULL)
        p4est_ghost_destroy(kinfo.forest.ghost)
        kinfo.forest.ghost = Ptr{p4est_ghost_t}(C_NULL)
    end
    p4est_destroy(p4est)
    kinfo.forest.p4est = Ptr{p4est_t}(C_NULL)
    p4est_connectivity_destroy(conn)
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
    conn = pointer(PointerWrapper(p4est).connectivity)
    if kinfo.forest.mesh != Ptr{p8est_mesh_t}(C_NULL)
        p8est_mesh_destroy(kinfo.forest.mesh)
        kinfo.forest.mesh = Ptr{p8est_mesh_t}(C_NULL)
    end
    if kinfo.forest.ghost != Ptr{p8est_ghost_t}(C_NULL)
        p8est_ghost_destroy(kinfo.forest.ghost)
        kinfo.forest.ghost = Ptr{p8est_ghost_t}(C_NULL)
    end
    p8est_destroy(p4est)
    kinfo.forest.p4est = Ptr{p8est_t}(C_NULL)
    p8est_connectivity_destroy(conn)
end

"""
$(TYPEDSIGNATURES)
"""
function finalize!(p4est::Ptr{p8est_t}, ka::KA)
    finalize_ghost!(ka.kdata.ghost.ghost_buffer)
    finalize_p4est!(p4est, ka)
    return nothing
end
