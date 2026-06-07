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
        dg*=0.5*(erf((c[i]+0.5*du[i])/vth)-erf((c[i]-0.5*du[i])/vth))
    end
    return dg
end
function maxwellian_refine_flag(midpoint,du,U,prim,kinfo::KInfo{DIM}) where{DIM}
    return maxwellian_dρ(midpoint,du,U,prim)> 1e-3*prim[1]
end

# ------------------------------------------------------------------------------------------
# Moment-weighted Löhner velocity-space refinement indicator (ADAPT_VS_MODE == :lohner).
#
# A model-free, single-grid interpolation-error indicator on the *distribution* `df`
# (not `cdf`): the normalized second velocity difference of the density and energy moments,
# evaluated against the same-or-coarser face neighbors (`vs_face_neighbor`), with the vacuum
# tail (`neighbor == 0`) treated as `f = 0`.  This is the velocity-space analogue of the
# physical-space Löhner sensor and replaces the magnitude/contribution refine flag, with
# `local_contribution_*` retained as a relative mass/energy floor.
# ------------------------------------------------------------------------------------------

"Löhner ripple-filter coefficient (relative noise floor)."
const VS_LOHNER_EPS_FLT = 1e-2
"Löhner absolute floor coefficient: structure below `EPS_ABS * scale` is ignored."
const VS_LOHNER_EPS_ABS = 1e-3
"Coarsen-eligibility fraction of the refine threshold (hysteresis)."
const VS_LOHNER_COARSEN_RATIO = 0.3

# Normalized non-uniform Löhner ratio from left/center/right values at center-to-neighbor
# distances `dsL`/`dsR`.  Same weighting as the physical-space `lohner_value`, so the second
# difference is correct across velocity cells of differing size (coarse/fine interfaces).
@inline function _lohner_ratio(L::Real, C::Real, R::Real, dsL::Real, dsR::Real, εflt::Real, εabs::Real, scale::Real)
    num = abs(dsR * L - (dsL + dsR) * C + dsL * R)
    den = dsR * abs(L - C) + dsL * abs(R - C) +
          εflt * (dsR * abs(L) + (dsL + dsR) * abs(C) + dsL * abs(R)) +
          εabs * scale * (dsL + dsR)
    return den > 0 ? num / den : 0.0
end

# Velocity-cell size along dimension `d`.
@inline _vs_cell_size(idx::VsNeighborIndex, vs::AbstractVsData, j::Integer, d::Integer) =
    @inbounds idx.h_fine[d] * (1 << (idx.maxlevel - vs.level[j]))

"""
$(TYPEDSIGNATURES)
Per-physical-cell reference magnitudes (peaks of the mass distribution `h` and, for
`NDF == 2`, the internal-energy distribution `b`) used as the absolute floor `scale` in the
Löhner ratio.
"""
function vs_lohner_scales(vs::AbstractVsData{DIM,NDF}) where {DIM,NDF}
    s1 = 0.0; s2 = 0.0
    @inbounds for j in 1:vs.vs_num
        a1 = abs(vs.df[j, 1]); a1 > s1 && (s1 = a1)
        if NDF == 2
            a2 = abs(vs.df[j, 2]); a2 > s2 && (s2 = a2)
        end
    end
    return s1, s2
end

"""
$(TYPEDSIGNATURES)
Moment-weighted Löhner indicator `η_i ∈ [0,1]` of velocity cell `i`: the max over each
velocity dimension of the normalized second difference of the distribution components
against the cell's face neighbors — the mass distribution `h = df[:,1]` and, for `NDF == 2`,
the independent internal-energy distribution `b = df[:,2]`.  The second difference is taken
on the raw distributions (no `c^2` weight, which would inject the spurious curvature of the
velocity weight itself).  `s1`/`s2` are the scales from [`vs_lohner_scales`](@ref).
"""
function vs_lohner_indicator(vs::AbstractVsData{DIM,NDF}, idx::VsNeighborIndex{DIM}, i::Integer, s1::Real, s2::Real) where {DIM,NDF}
    εflt = VS_LOHNER_EPS_FLT; εabs = VS_LOHNER_EPS_ABS
    @inbounds begin
        h0 = vs.df[i, 1]
        b0 = NDF == 2 ? vs.df[i, 2] : 0.0
        η = 0.0
        for d in 1:DIM
            hi = _vs_cell_size(idx, vs, i, d)
            Ln = vs_face_neighbor(idx, vs, i, d, -1)
            Rn = vs_face_neighbor(idx, vs, i, d, +1)
            # Center-to-center distance; the vacuum side (no neighbor) mirrors this cell.
            dsL = Ln == 0 ? hi : 0.5 * (hi + _vs_cell_size(idx, vs, Ln, d))
            dsR = Rn == 0 ? hi : 0.5 * (hi + _vs_cell_size(idx, vs, Rn, d))
            hL = Ln == 0 ? 0.0 : vs.df[Ln, 1]
            hR = Rn == 0 ? 0.0 : vs.df[Rn, 1]
            η = max(η, _lohner_ratio(hL, h0, hR, dsL, dsR, εflt, εabs, s1))
            if NDF == 2
                bL = Ln == 0 ? 0.0 : vs.df[Ln, 2]
                bR = Rn == 0 ? 0.0 : vs.df[Rn, 2]
                η = max(η, _lohner_ratio(bL, b0, bR, dsL, dsR, εflt, εabs, s2))
            end
        end
    end
    return η
end

# Carry the precomputed per-cell Löhner hint array alongside `df` through refine/coarsen so
# that `lnflag[index]` stays aligned with the mutating velocity grid (mirrors `cdf`).
function lohner_flag_refine_replace!(DIM::Integer, lnflag::AbstractVector{Bool}, index::Int)
    deleteat!(lnflag, index)
    for _ in 1:2^DIM
        insert!(lnflag, index, false)
    end
end
function lohner_flag_coarsen_replace!(DIM::Integer, lnflag::AbstractVector{Bool}, index::Int)
    for _ in 1:2^DIM-1
        deleteat!(lnflag, index)
    end
end