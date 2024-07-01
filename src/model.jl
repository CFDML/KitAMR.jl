ref_vhs_vis(Kn, alpha, omega) =
    5.0 * (alpha + 1.0) * (alpha + 2.0) * √π /
    (4.0 * alpha * (5.0 - 2.0 * omega) * (7.0 - 2.0 * omega)) * Kn
get_τ(prim::AV, μ::Real, ω::Real) = μ * 2.0 * prim[4]^(1 - ω) / prim[1]
