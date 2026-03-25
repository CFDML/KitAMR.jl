"""
$(TYPEDSIGNATURES)
"""
function adaptive_mesh_refinement!(p4est::P_pxest_t,ka::KA;ps_interval=40,vs_interval=80,partition_interval=40)
    ka.kinfo.status.residual.redundant_step>0&&(return nothing)
    res = maximum(ka.kinfo.status.residual.residual)
    converge_ratio = res/ka.kinfo.config.solver.TOLERANCE>100 ? 1 : Int(floor(100*ka.kinfo.config.solver.TOLERANCE/res))
    flag = false
    if ka.kinfo.config.solver.PS_DYNAMIC_AMR&&ka.kinfo.status.ps_adapt_step > ps_interval*converge_ratio
        ps_adaptive_mesh_refinement!(p4est,ka)
        flag = true;ka.kinfo.status.ps_adapt_step = 0
    end
    if ka.kinfo.config.solver.VS_DYNAMIC_AMR&&ka.kinfo.status.vs_adapt_step > vs_interval*converge_ratio
        vs_adaptive_mesh_refinement!(ka)
        flag = true;ka.kinfo.status.vs_adapt_step=0
    end
    if (ka.kinfo.config.solver.PS_DYNAMIC_AMR||ka.kinfo.config.solver.VS_DYNAMIC_AMR)&&ka.kinfo.status.partition_step>partition_interval*converge_ratio
        ps_partition!(p4est, ka)
        flag = true;ka.kinfo.status.partition_step = 0
    end
    if flag
        amr_recover!(p4est,ka)        
    end
    return nothing
end
function update_faces!(p4est::P_pxest_t, ka::KA)
    ka.kdata.field.faces = Vector{AbstractFace}(undef, 0)
    initialize_faces!(p4est, ka)
end

"""
$(TYPEDSIGNATURES)
Recover the ghost layers, neighbor relations, immersed boundaries, and faces after an AMR or partition process.
"""
function amr_recover!(p4est::P_pxest_t,ka::KA)
    update_ghost!(p4est, ka)
    update_neighbor!(p4est, ka)
    update_solid!(ka)
    update_faces!(p4est, ka)
    return nothing
end