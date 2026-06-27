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
function local_heatflux_contribution_ratio(w::AbstractVector, U::AbstractVector, midpoint::AbstractVector, df::AbstractVector, weight::Float64, ::KInfo{DIM,2}) where{DIM}
    c2 = 0.0
    @inbounds for d in 1:DIM
        c = midpoint[d] - U[d]
        c2 += c^2
    end
    _, _, heatflux = _relative_moment_ratios(
        w, U, c2, df[1] * weight, 0.5 * (c2 * df[1] + df[2]) * weight)
    return heatflux
end
function local_heatflux_contribution_ratio(w::AbstractVector, U::AbstractVector, midpoint::AbstractVector, df::AbstractVector, weight::Float64, ::KInfo{DIM,1}) where{DIM}
    c2 = 0.0
    @inbounds for d in 1:DIM
        c = midpoint[d] - U[d]
        c2 += c^2
    end
    _, _, heatflux = _relative_moment_ratios(
        w, U, c2, df[1] * weight, 0.5 * c2 * df[1] * weight)
    return heatflux
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
"Absolute floor coefficient used by velocity-space refinement indicators."
const VS_CRITERION_EPS_ABS = 1e-3

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
    return abs(dρ - Mmid * V) / max(dρ, VS_CRITERION_EPS_ABS * pref * V)
end

"""
$(TYPEDSIGNATURES)
Analytic initial-grid refine flag for velocity cell with center `midpoint`, size `du`.  Mirrors
the dynamic criterion using exact Maxwellian quantities: refine if the relative mass *or*
energy contribution exceeds [`MAXWELLIAN_INIT_FLOOR`](@ref), or if the analytic Maxwellian
quadrature error exceeds the internal `ADAPT_COEFFI_VS_INIT` default.  This initialization
criterion is independent of the dynamic velocity-space AMR pass.  `I0buf`/`I2buf` are reused
length-`DIM` scratch vectors.
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
    return maxwellian_quad_error(dρ, midpoint, du, U, prim, Val(DIM)) > ADAPT_COEFFI_VS_INIT
end

# Velocity-cell size along dimension `d`.
@inline _vs_cell_size(idx::VsNeighborIndex, vs::AbstractVsData, j::Integer, d::Integer) =
    @inbounds idx.h_fine[d] * (1 << (idx.maxlevel - vs.level[j]))

"""
$(TYPEDSIGNATURES)
Per-physical-cell reference magnitudes (peaks of the mass distribution `h` and, for
`NDF == 2`, the internal-energy distribution `b`) used as absolute floors in the LSR residual
normalization.
"""
function vs_lsr_scales(vs::AbstractVsData{DIM,NDF}) where {DIM,NDF}
    s1 = 0.0; s2 = 0.0
    @inbounds for j in 1:vs.vs_num
        a1 = abs(vs.df[j, 1]); a1 > s1 && (s1 = a1)
        if NDF == 2
            a2 = abs(vs.df[j, 2]); a2 > s2 && (s2 = a2)
        end
    end
    return s1, s2
end

# ------------------------------------------------------------------------------------------
# Local linear least-squares residual indicator.
#
# This sensor is dimension-unsplit: it fits a local affine model in the actual velocity
# coordinates around one velocity cell and measures the normalized residual.  Linear slopes
# are absorbed by the fit, so coarse/fine center offsets do not masquerade as directional
# second derivatives.  Boundary vacuum points are not injected; only existing face-neighbor
# 1-ring velocity cells participate in the local fit.
# ------------------------------------------------------------------------------------------

"Coarsen-eligibility fraction of the LSR refine threshold (hysteresis)."
const VS_LSR_COARSEN_RATIO = 0.3

@inline function _vs_root_ds(kinfo::KInfo{DIM}) where {DIM}
    return ntuple(d -> (kinfo.config.quadrature[2*d] - kinfo.config.quadrature[2*d - 1]) /
                      kinfo.config.vs_trees_num[d], DIM)
end

@inline _vs_local_lmax_floor(maxlevel::Integer) = maxlevel > 0 ? 1 : 0

function _maxwellian_haar_1d_max_rel(λ::Real, h_parent::Real)
    λ > 0 || return Inf
    h_parent > 0 || return 0.0
    σ = 1.0 / sqrt(2.0 * λ)
    a = 0.25 * h_parent
    xmax = max(8.0 * σ, 4.0 * h_parent)
    max_rel = 0.0
    # The density/prefactor cancels because the coefficient is normalized by the peak.
    @inbounds for i in 0:512
        x = xmax * i / 512
        rel = 0.5 * abs(exp(-λ * (x - a)^2) - exp(-λ * (x + a)^2))
        rel > max_rel && (max_rel = rel)
    end
    return max_rel
end

@inline vs_local_lmax_haar_threshold(::KInfo) = AMR_VS_HAAR_THRESHOLD

function analytic_maxwellian_local_lmax(prim::AbstractVector, kinfo::KInfo)
    maxlevel = kinfo.config.solver.AMR_VS_MAXLEVEL
    maxlevel <= 0 && return 0
    threshold = vs_local_lmax_haar_threshold(kinfo)
    floor_level = _vs_local_lmax_floor(maxlevel)
    ds0 = maximum(_vs_root_ds(kinfo))
    λ = prim[end]
    λ > 0 || return maxlevel
    for level in floor_level:maxlevel
        level == 0 && return 0
        h_parent = ds0 / 2.0^(level - 1)
        _maxwellian_haar_1d_max_rel(λ, h_parent) <= threshold && return level
    end
    return maxlevel
end

@inline function _vs_same_level_group(vs_data::AbstractVsData, first::Int, level::Int, nc::Int)
    first + nc - 1 <= vs_data.vs_num || return false
    @inbounds for g in 1:nc-1
        Int(vs_data.level[first + g]) == level || return false
    end
    return true
end

function _vs_sibling_groups(vs_data::AbstractVsData{DIM}, maxlevel::Integer) where {DIM}
    groups = Tuple{Int,Int}[]
    maxlevel <= 0 && return groups
    nc = 2^DIM
    align = zeros(Float64, maxlevel)
    index = 1
    @inbounds while index <= vs_data.vs_num
        level = Int(vs_data.level[index])
        if level > 0
            aligned = abs(mod(align[level], 1.0)) < 1.0e-12
            if aligned && _vs_same_level_group(vs_data, index, level, nc)
                push!(groups, (index, level))
                index += nc
                if level > 1
                    for l in 1:level-1
                        align[l] += 1 / 2^(DIM * (level - l))
                    end
                end
            else
                for l in 1:level
                    align[l] += 1 / 2^(DIM * (level - l + 1))
                end
                index += 1
            end
        else
            index += 1
        end
    end
    return groups
end

function _vs_component_scales(vs_data::AbstractVsData{DIM,NDF}) where {DIM,NDF}
    scales = Vector{Float64}(undef, NDF)
    @inbounds for k in 1:NDF
        s = 0.0
        for i in 1:vs_data.vs_num
            a = abs(vs_data.df[i, k])
            a > s && (s = a)
        end
        scales[k] = max(s, EPS)
    end
    return scales
end

@inline function _haar_child_sign(::Val{DIM}, child::Integer, mask::Integer) where {DIM}
    s = 1
    @inbounds for d in 1:DIM
        if (mask & (1 << (d - 1))) != 0
            s *= RMT[DIM][child][d]
        end
    end
    return s
end

function _vs_haar_group_rel_detail(vs_data::AbstractVsData{DIM,NDF}, first::Int,
                                   scales::AbstractVector) where {DIM,NDF}
    nc = 2^DIM
    max_rel = 0.0
    @inbounds for k in 1:NDF
        sumsq = 0.0
        for mask in 1:nc-1
            coeff = 0.0
            for child in 1:nc
                coeff += _haar_child_sign(Val(DIM), child, mask) *
                         vs_data.df[first + child - 1, k]
            end
            coeff /= nc
            sumsq += coeff * coeff
        end
        rel = sqrt(sumsq) / scales[k]
        rel > max_rel && (max_rel = rel)
    end
    return max_rel
end

@inline function _vs_heatflux_energy_density(vs_data::AbstractVsData{DIM,NDF},
                                             i::Integer,
                                             c2::Real) where {DIM,NDF}
    if NDF == 2
        return 0.5 * (c2 * vs_data.df[i, 1] + vs_data.df[i, 2])
    else
        return 0.5 * c2 * vs_data.df[i, 1]
    end
end

@inline function _vs_heatflux_integrand(vs_data::AbstractVsData{DIM,NDF},
                                        i::Integer,
                                        U::AbstractVector,
                                        dir::Integer) where {DIM,NDF}
    cdir = 0.0
    c2 = 0.0
    @inbounds for d in 1:DIM
        c = vs_data.midpoint[i, d] - U[d]
        d == dir && (cdir = c)
        c2 += c * c
    end
    return cdir * _vs_heatflux_energy_density(vs_data, i, c2)
end

function _vs_heatflux_scales(vs_data::AbstractVsData{DIM,NDF},
                             U::AbstractVector) where {DIM,NDF}
    scales = Vector{Float64}(undef, DIM)
    @inbounds for d in 1:DIM
        s = 0.0
        for i in 1:vs_data.vs_num
            a = abs(_vs_heatflux_integrand(vs_data, i, U, d))
            a > s && (s = a)
        end
        scales[d] = max(s, EPS)
    end
    return scales
end

function _vs_heatflux_haar_group_rel_detail(vs_data::AbstractVsData{DIM,NDF},
                                            first::Int,
                                            scales::AbstractVector,
                                            U::AbstractVector) where {DIM,NDF}
    nc = 2^DIM
    max_rel = 0.0
    @inbounds for d in 1:DIM
        sumsq = 0.0
        for mask in 1:nc-1
            coeff = 0.0
            for child in 1:nc
                coeff += _haar_child_sign(Val(DIM), child, mask) *
                         _vs_heatflux_integrand(vs_data, first + child - 1, U, d)
            end
            coeff /= nc
            sumsq += coeff * coeff
        end
        rel = sqrt(sumsq) / scales[d]
        rel > max_rel && (max_rel = rel)
    end
    return max_rel
end

function vs_haar_level_max_rel_detail(vs_data::AbstractVsData{DIM,NDF},
                                      level::Integer) where {DIM,NDF}
    level <= 0 && return 0.0
    walk_maxlevel = max(Int(level), maximum(Int.(vs_data.level)))
    scales = _vs_component_scales(vs_data)
    max_rel = 0.0
    for (first, group_level) in _vs_sibling_groups(vs_data, walk_maxlevel)
        group_level == level || continue
        rel = _vs_haar_group_rel_detail(vs_data, first, scales)
        rel > max_rel && (max_rel = rel)
    end
    return max_rel
end

function vs_heatflux_haar_level_max_rel_detail(vs_data::AbstractVsData{DIM,NDF},
                                               prim::AbstractVector,
                                               level::Integer) where {DIM,NDF}
    level <= 0 && return 0.0
    walk_maxlevel = max(Int(level), maximum(Int.(vs_data.level)))
    U = @view(prim[2:1+DIM])
    scales = _vs_heatflux_scales(vs_data, U)
    max_rel = 0.0
    for (first, group_level) in _vs_sibling_groups(vs_data, walk_maxlevel)
        group_level == level || continue
        rel = _vs_heatflux_haar_group_rel_detail(vs_data, first, scales, U)
        rel > max_rel && (max_rel = rel)
    end
    return max_rel
end

function vs_haar_virtual_coarsened_level_detail(vs_data::VsData{DIM,NDF}, ds,
                                                maxlevel::Integer,
                                                target_lmax::Integer) where {DIM,NDF}
    target_lmax <= 0 && return 0.0
    tmp = deepcopy(vs_data)
    coarsen_ok = [Int(level) > target_lmax for level in tmp.level]
    any(coarsen_ok) && coarsen_grid_stream!(tmp, coarsen_ok, ds, maxlevel)
    return vs_haar_level_max_rel_detail(tmp, target_lmax)
end

function vs_heatflux_haar_virtual_coarsened_level_detail(vs_data::VsData{DIM,NDF},
                                                        prim::AbstractVector,
                                                        ds,
                                                        maxlevel::Integer,
                                                        target_lmax::Integer) where {DIM,NDF}
    target_lmax <= 0 && return 0.0
    tmp = deepcopy(vs_data)
    coarsen_ok = [Int(level) > target_lmax for level in tmp.level]
    any(coarsen_ok) && coarsen_grid_stream!(tmp, coarsen_ok, ds, maxlevel)
    return vs_heatflux_haar_level_max_rel_detail(tmp, prim, target_lmax)
end

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
    den = sqrt(variation / wsum) + VS_CRITERION_EPS_ABS * max(abs(scale), abs(f0), EPS)
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
