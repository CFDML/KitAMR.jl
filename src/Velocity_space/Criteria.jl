"""
$(TYPEDSIGNATURES)
"""
function contribution_refine_flag(w::AbstractVector, U::AbstractVector, midpoint::AbstractVector, df::AbstractVector, weight::Float64, vr::Velocity_Resolution, kinfo::KInfo)
    return local_contribution_refine_flag(w,U,midpoint,df,weight,kinfo)||global_contribution_refine_flag(U,midpoint,df,weight,vr,kinfo)
end
"""
$(TYPEDSIGNATURES)
"""
function contribution_coarsen_flag(w::AbstractVector, U::AbstractVector, midpoint::AbstractMatrix, df::AbstractMatrix, weight::AbstractVector, vr::Velocity_Resolution, kinfo::KInfo)
    return local_contribution_coarsen_flag(w,U,midpoint,df,weight,kinfo)&&global_contribution_coarsen_flag(U,midpoint,df,weight,vr,kinfo)
end





function local_contribution_refine_flag(w::AbstractVector, U::AbstractVector, midpoint::AbstractVector, df::AbstractVector, weight::Float64, kinfo::KInfo{DIM,2}) where{DIM}
    max(abs(0.5 * (sum((U - midpoint) .^ 2) * df[1] + df[2])* weight) /
    (w[end]  - 0.5 * w[1] * sum((U) .^ 2)),df[1]*weight/w[1]) > kinfo.config.solver.ADAPT_COEFFI_VS_LOCAL ? true : false
end
function local_contribution_refine_flag(w::AbstractVector, U::AbstractVector, midpoint::AbstractVector, df::AbstractVector, weight::Float64, kinfo::KInfo{DIM,1}) where{DIM}
    max(abs(0.5 * sum((U - midpoint) .^ 2) * df[1]* weight) /
    (w[end] - 0.5 * w[1] * sum((U) .^ 2)),df[1]*weight/w[1]) > kinfo.config.solver.ADAPT_COEFFI_VS_LOCAL ? true : false
end
function local_contribution_coarsen_flag(w::AbstractVector, U::AbstractVector, midpoint::AbstractMatrix, df::AbstractMatrix, weight::AbstractVector, kinfo::KInfo{DIM,2}) where{DIM}
    for i = 1:2^DIM
        !local_contribution_coarsen_flag(w, U, @view(midpoint[i, :]), @view(df[i, :]), weight[i], kinfo) &&
            return false
    end
    return true
end
function local_contribution_coarsen_flag(w::AbstractVector, U::AbstractVector, midpoint::AbstractMatrix, df::AbstractMatrix, weight::AbstractVector, kinfo::KInfo{DIM,1}) where{DIM}
    for i = 1:2^DIM
        !local_contribution_coarsen_flag(w, U, @view(midpoint[i, :]), @view(df[i, :]), weight[i], kinfo) &&
            return false
    end
    return true
end
function local_contribution_coarsen_flag(w::AbstractVector, U::AbstractVector, midpoint::AbstractVector, df::AbstractVector, weight::Float64, kinfo::KInfo{DIM,2}) where{DIM}
    max(0.5 * (sum((U - midpoint) .^ 2) * df[1] + df[2])* weight /
    (w[end] - 0.5 * w[1] * sum((U) .^ 2)),df[1]*weight/w[1]) < kinfo.config.solver.ADAPT_COEFFI_VS_LOCAL/ 2^DIM
end
function local_contribution_coarsen_flag(w::AbstractVector, U::AbstractVector, midpoint::AbstractVector, df::AbstractVector, weight::Float64, kinfo::KInfo{DIM,1}) where{DIM}
    max(0.5 * (sum((U - midpoint) .^ 2) * df[1])* weight /
    (w[end] - 0.5 * w[1] * sum((U) .^ 2)),df[1]*weight/w[1]) < kinfo.config.solver.ADAPT_COEFFI_VS_LOCAL/ 2^DIM
end





function global_contribution_refine_flag(U::AbstractVector, midpoint::AbstractVector, df::AbstractVector, weight::Float64, vr::Velocity_Resolution, kinfo::KInfo{DIM,2}) where{DIM}
    df[1]*weight > kinfo.config.solver.ADAPT_COEFFI_VS_GLOBAL*vr.density || 0.5 * (sum((U - midpoint) .^ 2) * df[1] + df[2])* weight>vr.energy*kinfo.config.solver.ADAPT_COEFFI_VS_GLOBAL
end
function global_contribution_refine_flag(U::AbstractVector, midpoint::AbstractVector, df::AbstractVector, weight::Float64, vr::Velocity_Resolution, kinfo::KInfo{DIM,1}) where{DIM}
    df[1]*weight > vr.density*kinfo.config.solver.ADAPT_COEFFI_VS_GLOBAL || 0.5 * (sum((U - midpoint) .^ 2) * df[1])* weight>vr.energy*kinfo.config.solver.ADAPT_COEFFI_VS_GLOBAL
end

function global_contribution_coarsen_flag(U::AbstractVector, midpoint::AbstractMatrix, df::AbstractMatrix, weight::AbstractVector, vr::Velocity_Resolution, kinfo::KInfo{DIM,2}) where{DIM}
    for i = 1:2^DIM
        !global_contribution_coarsen_flag(U, @view(midpoint[i, :]), @view(df[i, :]), weight[i], vr, kinfo) &&
            return false
    end
    return true
end
function global_contribution_coarsen_flag(U::AbstractVector, midpoint::AbstractMatrix, df::AbstractMatrix, weight::AbstractVector, vr::Velocity_Resolution, kinfo::KInfo{DIM,1}) where{DIM}
    for i = 1:2^DIM
        !global_contribution_coarsen_flag(U, @view(midpoint[i, :]), @view(df[i, :]), weight[i], vr, kinfo) &&
            return false
    end
    return true
end
function global_contribution_coarsen_flag(U::AbstractVector, midpoint::AbstractVector, df::AbstractVector, weight::Float64, vr::Velocity_Resolution, kinfo::KInfo{DIM,2}) where{DIM}
    df[1]*weight < kinfo.config.solver.ADAPT_COEFFI_VS_GLOBAL*vr.density/2^DIM && 0.5 * (sum((U - midpoint) .^ 2) * df[1] + df[2])* weight<vr.energy*kinfo.config.solver.ADAPT_COEFFI_VS_GLOBAL/2^DIM
end
function global_contribution_coarsen_flag(U::AbstractVector, midpoint::AbstractVector, df::AbstractVector, weight::Float64, vr::Velocity_Resolution, kinfo::KInfo{DIM,1}) where{DIM}
    df[1]*weight < vr.density*kinfo.config.solver.ADAPT_COEFFI_VS_GLOBAL/2^DIM && 0.5 * (sum((U - midpoint) .^ 2) * df[1])* weight<vr.energy*kinfo.config.solver.ADAPT_COEFFI_VS_GLOBAL/2^DIM
end

function macro_estimate_refine_flag(prim::AbstractVector,U,midpoint,ds,level)
    return false
    if norm(midpoint-U)-norm(ds)/2^(level+1)<√(1/prim[end])
        return true
    else
        return false
    end
end
function macro_estimate_IB_refine_flag(prim,bc,U,midpoint,ds,level)
    Ub = bc[2:length(midpoint)+1]
    if (norm(midpoint-U)-norm(ds)/2^(level+1)<√(1/prim[end])||norm(midpoint-Ub)-norm(ds)/2^(level+1)<√(1/bc[end]))
        return true
    else
        return false
    end
end

"""
$(SIGNATURES)
Analytical integral of the Maxwellian distribution over a Cartesian grid.
"""
function maxwellian_dρ(midpoint,du,U,prim)
    c = midpoint-U
    dg = prim[1]
    vth = 1.0/√prim[end]
    for i in eachindex(midpoint)
        dg*=0.5/vth*(erf((c[i]+0.5*du[i])/vth)-erf((c[i]-0.5*du[i])/vth))
    end
    return dg
end
function maxwellian_refine_flag(midpoint,du,U,prim,kinfo::KInfo{DIM}) where{DIM}
    return maxwellian_dρ(midpoint,du,U,prim)> 1e-3*prim[1]
end