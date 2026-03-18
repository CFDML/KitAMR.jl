"""
$(TYPEDSIGNATURES)
"""
function adaptive_mesh_refinement!(ps4est::P_pxest_t,amr::KitAMR_Data;ps_interval=40,vs_interval=80,partition_interval=40)
    amr.global_data.status.residual.redundant_step>0&&(return nothing)
    res = maximum(amr.global_data.status.residual.residual)
    converge_ratio = res/TOLERANCE>100 ? 1 : Int(floor(100*TOLERANCE/res))
    flag = false
    if amr.global_data.config.solver.PS_DYNAMIC_AMR&&amr.global_data.status.ps_adapt_step > ps_interval*converge_ratio
        ps_adaptive_mesh_refinement!(ps4est,amr)
        flag = true;amr.global_data.status.ps_adapt_step = 0
    end
    if amr.global_data.config.solver.VS_DYNAMIC_AMR&&amr.global_data.status.vs_adapt_step > vs_interval*converge_ratio
        vs_adaptive_mesh_refinement!(amr)
        flag = true;amr.global_data.status.vs_adapt_step=0
    end
    if (amr.global_data.config.solver.PS_DYNAMIC_AMR||amr.global_data.config.solver.VS_DYNAMIC_AMR)&&amr.global_data.status.partition_step>partition_interval*converge_ratio
        ps_partition!(ps4est, amr)
        flag = true;amr.global_data.status.partition_step = 0
    end
    if flag
        amr_recover!(ps4est,amr)        
    end
    return nothing
end
function update_faces!(p4est::P_pxest_t, amr::KitAMR_Data)
    amr.field.faces = Vector{AbstractFace}(undef, 0)
    initialize_faces!(p4est, amr)
end

"""
$(TYPEDSIGNATURES)
Recover the ghost layers, neighbor relations, immersed boundaries, and faces after an AMR or partition process.
"""
function amr_recover!(ps4est::P_pxest_t,amr::KitAMR_Data)
    update_ghost!(ps4est, amr)
    update_neighbor!(ps4est, amr)
    update_solid!(amr)
    update_faces!(ps4est, amr)
    return nothing
end