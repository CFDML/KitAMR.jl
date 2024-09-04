function get_conserved_3D(prim::AbstractVector, γ::Real)
    w = Vector{Float64}(undef, 5)
    @inbounds w[1] = prim[1]
    @inbounds w[2] = prim[1] * prim[2]
    @inbounds w[3] = prim[1] * prim[3]
    @inbounds w[4] = prim[1] * prim[4]
    @inbounds w[5] =
        0.5 * prim[1] / prim[5] / (γ - 1.0) +
        0.5 * prim[1] * (prim[2]^2 + prim[3]^2 + prim[4]^2)
    return w
end
function get_prim_3D(w::AbstractVector, γ::Real)
    prim = Vector{Float64}(undef, 5)
    @inbounds prim[1] = w[1]
    @inbounds prim[2] = w[2] / w[1]
    @inbounds prim[3] = w[3] / w[1]
    @inbounds prim[4] = w[4] / w[1]
    @inbounds prim[5] = 0.5 * w[1] / (γ - 1.0) / (w[5] - 0.5 * (w[2]^2 + w[3]^2 + w[4]^2) / w[1])
    return prim
end
function micro_slope_3D(sw::AbstractVector{T}, prim::AbstractVector, K::Real) where {T}
    a = Vector{T}(undef, 5)
    @inbounds a[5] =
        4.0 * prim[5]^2 / (K + 3.0) / prim[1] * (
            2.0 * sw[5] - 2.0 * prim[2] * sw[2] - 2.0 * prim[3] * sw[3] -
            2.0 * prim[4] * sw[4] +
            sw[1] * (prim[2]^2 + prim[3]^2 + prim[4]^2 - 0.5 * (K + 3.0) / prim[5])
        )
    @inbounds a[4] = 2 * prim[5] / prim[1] * (sw[4] - prim[4] * sw[1]) - prim[4] * a[5]
    @inbounds a[3] = 2 * prim[5] / prim[1] * (sw[3] - prim[3] * sw[1]) - prim[3] * a[5]
    @inbounds a[2] = 2 * prim[5] / prim[1] * (sw[2] - prim[2] * sw[1]) - prim[2] * a[5]
    @inbounds a[1] =
        sw[1] / prim[1] - prim[2] * a[2] - prim[3] * a[3] - prim[4] * a[4] -
        0.5 * (prim[2]^2 + prim[3]^2 + prim[4]^2 + 0.5 * (K + 3.0) / prim[5]) * a[5]
    return a
end