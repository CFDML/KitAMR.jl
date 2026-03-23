"""
$(TYPEDSIGNATURES)
Reference viscousity obtained by VHS molecules model.
"""
function ref_vhs_vis(Kn, omega, alpha = 1.0, T_ref = 1.0)
    μ_0 =
        5.0 * (alpha + 1.0) * (alpha + 2.0) * √π /
        (4.0 * alpha * (5.0 - 2.0 * omega) * (7.0 - 2.0 * omega)) * Kn
    return μ_0*(T_ref)^(0.5-omega)
end
"""
$(TYPEDSIGNATURES)
Calculate collisional time.
"""
get_τ(prim::AbstractVector, μ::Real, ω::Real) = μ * 2.0 * prim[end]^(1 - ω) / prim[1]
