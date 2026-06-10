"""
$(TYPEDSIGNATURES)
"""
function adaptive_mesh_refinement!(
    p4est::P_pxest_t,
    ka::KA;
    ps_interval = 40,
    vs_interval = 80,
    partition_interval = 40,
    ps_recursive = false,
    vs_balance = false,
)
    ka.kinfo.status.residual.redundant_step>0&&(return nothing)
    res = maximum(ka.kinfo.status.residual.residual)
    converge_ratio =
        res/ka.kinfo.config.solver.TOLERANCE>100 ? 1 :
        Int(floor(100*ka.kinfo.config.solver.TOLERANCE/res))
    ps_changed = false
    vs_changed = false
    partition_changed = false
    if ka.kinfo.config.solver.PS_DYNAMIC_AMR&&ka.kinfo.status.ps_adapt_step >
                                              ps_interval*converge_ratio
        ps_adaptive_mesh_refinement!(p4est, ka; recursive = ps_recursive)
        ps_changed = true;
        ka.kinfo.status.ps_adapt_step = 1
    end
    if ka.kinfo.config.solver.VS_DYNAMIC_AMR&&ka.kinfo.status.vs_adapt_step >
                                              vs_interval*converge_ratio
        if vs_balance && ps_changed
            update_ghost!(p4est, ka)
            update_neighbor!(p4est, ka)
        end
        vs_changed = vs_adaptive_mesh_refinement!(ka; vs_balance = vs_balance)
        ka.kinfo.status.vs_adapt_step=1
    end
    if (ka.kinfo.config.solver.PS_DYNAMIC_AMR||ka.kinfo.config.solver.VS_DYNAMIC_AMR)&&ka.kinfo.status.partition_step>partition_interval*converge_ratio&&(
        ps_changed||vs_changed
    )
        ps_partition!(p4est, ka)
        partition_changed = true;
        ka.kinfo.status.partition_step = 1
    end
    if ps_changed || partition_changed
        amr_recover!(p4est, ka; topology_changed = true, velocity_changed = vs_changed)
    elseif vs_changed
        amr_recover!(p4est, ka; topology_changed = false, velocity_changed = true)
    end
    return nothing
end
function update_faces!(p4est::P_pxest_t, ka::KA)
    empty!(ka.kdata.field.faces)
    initialize_faces!(p4est, ka)
end

_has_immersed_boundaries(ka::KA) = !isempty(ka.kinfo.config.IB)

function _finish_immersed_boundary_recover!(ka::KA)
    if _has_immersed_boundaries(ka)
        initialize_immersed_boundaries!(ka)
    else
        empty!(ka.kdata.field.immersed_boundaries)
    end
    return nothing
end

"""
$(TYPEDSIGNATURES)
Recover the ghost layers, neighbor relations, immersed boundaries, and faces after an AMR or partition process.
"""
function amr_recover!(
    p4est::P_pxest_t,
    ka::KA;
    topology_changed::Bool = true,
    velocity_changed::Bool = true,
)
    if topology_changed
        update_ghost!(p4est, ka)
        update_neighbor!(p4est, ka)
    elseif velocity_changed && MPI.Comm_size(MPI.COMM_WORLD) > 1
        vs_ghost_exchange!(p4est, ka)
    end

    has_ib = _has_immersed_boundaries(ka)
    has_ib && update_solid!(ka)

    if topology_changed || has_ib || (velocity_changed && MPI.Comm_size(MPI.COMM_WORLD) > 1)
        update_faces!(p4est, ka)
    end
    _finish_immersed_boundary_recover!(ka)
    return nothing
end
