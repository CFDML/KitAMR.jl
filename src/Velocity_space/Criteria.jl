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

@inline function _local_moment_scales(w::AbstractVector, U::AbstractVector)
    u2 = 0.0
    @inbounds for d in eachindex(U)
        u2 += U[d]^2
    end
    rho = abs(w[1]) + EPS
    e_int = abs(w[end] - 0.5 * w[1] * u2) + EPS
    qscale = e_int * sqrt(max(e_int / rho, EPS)) + EPS
    return rho, e_int, qscale
end

@inline function _relative_moment_ratios(w::AbstractVector, U::AbstractVector, c2::Real,
                                         mass_cell::Real, energy_cell::Real)
    rho, e_int, qscale = _local_moment_scales(w, U)
    mass = abs(mass_cell) / rho
    energy = abs(energy_cell) / e_int
    heatflux = sqrt(c2) * abs(energy_cell) / qscale
    return mass, energy, heatflux
end

function local_contribution_ratio(w::AbstractVector, U::AbstractVector, midpoint::AbstractVector, df::AbstractVector, weight::Float64, ::KInfo{DIM,2}) where{DIM}
    c2 = 0.0
    @inbounds for d in 1:DIM
        c = midpoint[d] - U[d]
        c2 += c^2
    end
    mass, energy, heatflux = _relative_moment_ratios(
        w, U, c2, df[1] * weight, 0.5 * (c2 * df[1] + df[2]) * weight)
    return max(mass, energy, heatflux)
end
function local_contribution_ratio(w::AbstractVector, U::AbstractVector, midpoint::AbstractVector, df::AbstractVector, weight::Float64, ::KInfo{DIM,1}) where{DIM}
    c2 = 0.0
    @inbounds for d in 1:DIM
        c = midpoint[d] - U[d]
        c2 += c^2
    end
    mass, energy, heatflux = _relative_moment_ratios(
        w, U, c2, df[1] * weight, 0.5 * c2 * df[1] * weight)
    return max(mass, energy, heatflux)
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
# ------------------------------------------------------------------------------------------
# Analytic initial-grid refinement criterion.
#
# At initialization the distribution is exactly Maxwellian, so the cell curvature and the
# relative mass/energy contribution are evaluated from the *analytic* Maxwellian (exact
# erf-based cell integrals + analytic neighbor values) rather than from discrete samples.
# This is robust to a coarse starting grid: a cell that overlaps an important region (peak or
# rapidly-varying tail) is detected and refined even when no velocity node yet lies there,
# and recursive refinement then homes in.
# ------------------------------------------------------------------------------------------

"Relative mass/energy floor for the analytic initial-grid criterion."
const MAXWELLIAN_INIT_FLOOR = 1e-3

# 1-D integrals over [c-du/2, c+du/2] of exp(-λs²) and s²exp(-λs²) (analytic, erf-based).
@inline function _maxwellian_1d(c::Real, du::Real, λ::Real)
    sqλ = sqrt(λ); a = c - 0.5 * du; b = c + 0.5 * du
    I0 = 0.5 * sqrt(π / λ) * (erf(sqλ * b) - erf(sqλ * a))
    I2 = (a * exp(-λ * a^2) - b * exp(-λ * b^2)) / (2λ) + I0 / (2λ)
    return I0, I2
end

"""
$(TYPEDSIGNATURES)
Analytic under-resolution indicator of the Maxwellian over a cell of size `du`: the *relative*
exact-minus-midpoint quadrature error of the cell mass, `|∫M_h − M_h(mid)·V| / ∫M_h`, where
`dρ = ∫M_h` is the exact cell integral.  Because it compares the *exact* integral (which sees
the peak wherever it lies inside the cell) against the single midpoint sample, it cannot miss
sub-cell structure the way a `±du` finite difference can: it is O(1) once the cell is coarser
than the thermal width — regardless of how the cell straddles the peak — and is normalized by
the cell's own mass so a coarse cell carrying significant mass is flagged rather than diluted
by the (large) peak-times-volume.  The denominator carries a deep-tail floor so cells with
negligible mass (where both terms vanish) read ≈0 instead of `0/0`.
"""
function maxwellian_quad_error(dρ, midpoint, du, U, prim, ::Val{DIM}) where {DIM}
    λ = prim[end]; ρ = prim[1]
    pref = ρ * (λ / π)^(DIM / 2)
    V = 1.0; sumc2 = 0.0
    @inbounds for d in 1:DIM
        V *= du[d]
        sumc2 += (midpoint[d] - U[d])^2
    end
    Mmid = pref * exp(-λ * sumc2)
    return abs(dρ - Mmid * V) / max(dρ, VS_LOHNER_EPS_ABS * pref * V)
end

"""
$(TYPEDSIGNATURES)
Analytic initial-grid refine flag for velocity cell with center `midpoint`, size `du`.  Mirrors
the dynamic criterion using exact Maxwellian quantities: refine if the relative mass *or*
energy contribution exceeds [`MAXWELLIAN_INIT_FLOOR`](@ref), or if the analytic Maxwellian
quadrature error exceeds `ADAPT_COEFFI_VS_INIT`.  This initialization criterion is independent
of the dynamic velocity-space AMR mode.  `I0buf`/`I2buf` are reused length-`DIM` scratch vectors.
"""
function maxwellian_refine_flag(midpoint, du, U, prim, kinfo::KInfo{DIM,NDF}, I0buf, I2buf) where {DIM,NDF}
    λ = prim[end]; ρ = prim[1]
    @inbounds for d in 1:DIM
        I0buf[d], I2buf[d] = _maxwellian_1d(midpoint[d] - U[d], du[d], λ)
    end
    pref = ρ * (λ / π)^(DIM / 2)
    prodI0 = 1.0
    @inbounds for d in 1:DIM
        prodI0 *= I0buf[d]
    end
    dρ = pref * prodI0                            # exact cell mass
    dc2 = 0.0                                      # exact ∫ c² M_h over the cell
    @inbounds for d in 1:DIM
        p = 1.0
        for e in 1:DIM
            e != d && (p *= I0buf[e])
        end
        dc2 += I2buf[d] * p
    end
    dc2 *= pref
    if NDF == 2
        K = kinfo.config.gas.K
        dE = 0.5 * (dc2 + (K / (2λ)) * dρ)
        Eint = ρ * (DIM + K) / (4λ)
    else
        dE = 0.5 * dc2
        Eint = ρ * DIM / (4λ)
    end
    max(dρ / ρ, dE / Eint) > MAXWELLIAN_INIT_FLOOR && return true
    return maxwellian_quad_error(dρ, midpoint, du, U, prim, Val(DIM)) > kinfo.config.solver.ADAPT_COEFFI_VS_INIT
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

# ------------------------------------------------------------------------------------------
# Local linear least-squares residual indicator (ADAPT_VS_MODE == :lsr).
#
# This sensor is dimension-unsplit: it fits a local affine model in the actual velocity
# coordinates around one velocity cell and measures the normalized residual.  Linear slopes
# are absorbed by the fit, so coarse/fine center offsets do not masquerade as directional
# second derivatives.  Boundary vacuum points are not injected; only existing face-neighbor
# 1-ring velocity cells participate in the local fit.
# ------------------------------------------------------------------------------------------

@inline function _push_lsr_neighbor!(buf::Vector{Int}, nb::Integer, center::Integer)
    nb == 0 && return nothing
    nb == center && return nothing
    @inbounds for x in buf
        x == nb && return nothing
    end
    push!(buf, Int(nb))
    return nothing
end

function _lsr_neighbors!(buf::Vector{Int}, vs::AbstractVsData{DIM}, idx::VsNeighborIndex{DIM}, i::Integer) where {DIM}
    empty!(buf)
    @inbounds for d in 1:DIM
        _push_lsr_neighbor!(buf, vs_face_neighbor(idx, vs, i, d, -1), i)
        _push_lsr_neighbor!(buf, vs_face_neighbor(idx, vs, i, d, 1), i)
    end
    return buf
end

function _lsr_component_indicator(
    vs::AbstractVsData{DIM,NDF},
    idx::VsNeighborIndex{DIM},
    i::Integer,
    k::Integer,
    scale::Real,
    neighbors::Vector{Int},
    normal::AbstractMatrix{Float64},
    rhs::AbstractVector{Float64},
    x::AbstractVector{Float64},
) where {DIM,NDF}
    length(neighbors) < DIM + 1 && return 0.0
    fill!(normal, 0.0)
    fill!(rhs, 0.0)

    f0 = vs.df[i, k]
    wsum = 0.0
    variation = 0.0
    @inbounds for nb in neighbors
        dist2 = 0.0
        for d in 1:DIM
            h = _vs_cell_size(idx, vs, i, d)
            x[d] = (vs.midpoint[nb, d] - vs.midpoint[i, d]) / h
            dist2 += x[d]^2
        end
        dist2 <= EPS && continue
        wt = 1.0 / (1.0 + dist2)
        df = vs.df[nb, k] - f0
        wsum += wt
        variation += wt * df^2
        for a in 1:DIM
            rhs[a] += wt * x[a] * df
            for b in 1:DIM
                normal[a, b] += wt * x[a] * x[b]
            end
        end
    end
    wsum <= 0.0 && return 0.0

    traceN = 0.0
    @inbounds for d in 1:DIM
        traceN += normal[d, d]
    end
    ridge = 1e-12 * max(traceN, 1.0)
    @inbounds for d in 1:DIM
        normal[d, d] += ridge
    end

    gradient = try
        Symmetric(normal) \ rhs
    catch
        return 0.0
    end

    residual = 0.0
    @inbounds for nb in neighbors
        dist2 = 0.0
        for d in 1:DIM
            h = _vs_cell_size(idx, vs, i, d)
            x[d] = (vs.midpoint[nb, d] - vs.midpoint[i, d]) / h
            dist2 += x[d]^2
        end
        dist2 <= EPS && continue
        wt = 1.0 / (1.0 + dist2)
        pred = 0.0
        for d in 1:DIM
            pred += gradient[d] * x[d]
        end
        r = (vs.df[nb, k] - f0) - pred
        residual += wt * r^2
    end

    num = sqrt(residual / wsum)
    den = sqrt(variation / wsum) + VS_LOHNER_EPS_ABS * max(abs(scale), abs(f0), EPS)
    return min(num / den, 1.0)
end

function vs_lsr_indicator(
    vs::AbstractVsData{DIM,NDF},
    idx::VsNeighborIndex{DIM},
    i::Integer,
    s1::Real,
    s2::Real,
    neighbors::Vector{Int},
    normal::AbstractMatrix{Float64},
    rhs::AbstractVector{Float64},
    x::AbstractVector{Float64},
) where {DIM,NDF}
    _lsr_neighbors!(neighbors, vs, idx, i)
    η = _lsr_component_indicator(vs, idx, i, 1, s1, neighbors, normal, rhs, x)
    if NDF == 2
        η = max(η, _lsr_component_indicator(vs, idx, i, 2, s2, neighbors, normal, rhs, x))
    end
    return η
end
