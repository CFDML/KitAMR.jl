function discrete_maxwell_3D1F(
    u::AbstractVector{T},
    v::AbstractVector,
    w::AbstractVector,
    prim::AbstractVector,
) where {T}
    M = Vector{T}(undef, length(u))
    @inbounds @. M =
        prim[1] *
        (prim[5] / π)^(3 / 2) *
        exp(-prim[5] * ((u - prim[2])^2 + (v - prim[3])^2 + (w - prim[4])^2))
    return M
end
function shakhov_part_3D1F(
    u::AbstractVector{T},
    v::AbstractVector,
    w::AbstractVector,
    M::AbstractVector,
    prim::AbstractVector,
    qf::AbstractVector,
    Pr::Real,
) where {T}
    M⁺ = Vector{T}(undef, length(u))
    @inbounds @. M⁺ =
        0.8 * (1 - Pr) * prim[5]^2 / prim[1] *
        ((u - prim[2]) * qf[1] + (v - prim[3]) * qf[2] + (w - prim[4]) * qf[3]) *
        (2 * prim[5] * ((u - prim[2])^2 + (v - prim[3])^2 + (w - prim[4])^2) - 5) *
        M
    return M⁺
end
function heat_flux_3D1F(
    u::AbstractVector{T},
    v::AbstractVector,
    w::AbstractVector,
    h::AbstractVector,
    prim::AbstractVector,
    weight::AbstractVector,
) where{T}
    q = Vector{T}(undef, 3)
    @inbounds q[1] =
        0.5 * (sum(
            @. weight *
               (u - prim[2]) *
               ((u - prim[2])^2 + (v - prim[3])^2 + (w - prim[4])^2) *
               h
        ))
    @inbounds q[2] =
        0.5 * (sum(
            @. weight *
               (v - prim[3]) *
               ((u - prim[2])^2 + (v - prim[3])^2 + (w - prim[4])^2) *
               h
        ))
    @inbounds q[3] =
        0.5 * (sum(
            @. weight *
               (w - prim[4]) *
               ((u - prim[2])^2 + (v - prim[3])^2 + (w - prim[4])^2) *
               h
        ))
    return q
end
function maxwellian_density_3D1F(
    u::AbstractVector,
    v::AbstractVector,
    w::AbstractVector,
    h::AbstractVector,
    prim::AbstractVector,
    weight::AbstractVector,
    Θ::AbstractVector,
    vn::AbstractVector,
)
    @inbounds SF = sum(@. weight * vn * h * (1.0 - Θ))
    @inbounds SG =
        (prim[5] / π)^(3 / 2) * sum(
            @. weight *
               vn *
               exp(-prim[5] * ((u - prim[2])^2 + (v - prim[3])^2 + (w - prim[4])^2)) *
               Θ
        )
    return -SF / SG
end
function macro_flux_3D1F(
    u::AbstractVector{T},
    v::AbstractVector,
    w::AbstractVector,
    h::AbstractVector,
    weight::AbstractVector,
    vn::AbstractVector,
) where {T}
    flux = Vector{T}(undef, 5)
    @inbounds flux[1] = sum(@. weight * vn * h)
    @inbounds flux[2] = sum(@. weight * vn * u * h)
    @inbounds flux[3] = sum(@. weight * vn * v * h)
    @inbounds flux[4] = sum(@. weight * vn * w * h)
    @inbounds flux[5] = 0.5 * sum(@. weight * vn * (u^2 + v^2 + w^2) * h)
    return flux
end
function micro_to_macro_3D1F(
    u::AbstractVector{T},
    v::AbstractVector,
    w::AbstractVector,
    f::AbstractVector,
    weight::AbstractVector,
) where {T}
    W = Vector{T}(undef, 5)
    @inbounds W[1] = sum(@. weight * f)
    @inbounds W[2] = sum(@. weight * u * f)
    @inbounds W[3] = sum(@. weight * v * f)
    @inbounds W[4] = sum(@. weight * w * f)
    @inbounds W[5] = 0.5 * sum(@. weight * (u^2 + v^2 + w^2) * f)
    return W
end

function moment_u_3D1F(U::T,V::Real,W::Real,λ::Real,n::Integer,m::Integer,δ::Integer) where{T} # n>=1 calculate <u^n> <v^m> <w^δ> <u^n>_>0 <u^n>_<0 
    Mu_L = Vector{T}(undef, n + 1)
    Mu = Vector{T}(undef, n + 1)
    Mu_R = Vector{T}(undef, n + 1)
    @inbounds begin
        Mu_L[1] = 0.5 * erfc(-√(λ) * U)
        Mu_L[2] =
            U * Mu_L[1] -  0.5 * exp(-λ * U^2) / (√(π * λ))
        Mu_R[1] = 0.5 * erfc(√(λ) * U)
        Mu_R[2] =
            U * Mu_R[1] + 0.5 * exp(-λ * U^2) / (√(π * λ))
        for i = 1:n-1
            Mu_L[i+2] = U * Mu_L[i+1] + 0.5 * i / λ * Mu_L[i]
            Mu_R[i+2] = U * Mu_R[i+1] + 0.5 * i / λ * Mu_R[i]
        end
        Mu = @. Mu_L + Mu_R
        Mv = Vector{T}(undef, m + 1)
        Mv[1] = one(T)
        Mv[2] = V
        for i = 1:m-1
            Mv[i+2] = V * Mv[i+1] + 0.5 * i * Mv[i] / λ
        end
        Mw = Vector{T}(undef, δ + 1)
        Mw[1] = one(T)
        Mw[2] = W
        for i = 1:δ-1
            Mw[i+2] = W * Mw[i+1] + 0.5 * i * Mw[i] / λ
        end
    end
    return (Mu, Mv, Mw, Mu_L, Mu_R)
end

function moment_uv_3D1F(Mu::AbstractVector{T},Mv,Mw,α::Integer, β::Integer, δ::Integer) where{T}# calc <u^αv^βw^δψ>
    Muv = Vector{T}(undef, 5)
    @inbounds begin
        Muv[1] = Mu[α+1] * Mv[β+1] * Mw[δ+1]
        Muv[2] = Mu[α+2] * Mv[β+1] * Mw[δ+1]
        Muv[3] = Mu[α+1] * Mv[β+2] * Mw[δ+1]
        Muv[4] = Mu[α+1] * Mv[β+1] * Mw[δ+2]
        Muv[5] =
            0.5 * (
                Mu[α+3] * Mv[β+1] * Mw[δ+1] +
                Mu[α+1] * Mv[β+3] * Mw[δ+1] +
                Mu[α+1] * Mv[β+1] * Mw[δ+3]
            )
    end
    return Muv
end

function moment_au_3D1F(a, Mu, Mv, Mw, α::Integer, β::Integer, δ::Integer)# calc <au^αv^βw^δψ> 
    @inbounds a[1] * moment_uv_3D1F(Mu, Mv, Mw, α, β, δ) +
              a[2] * moment_uv_3D1F(Mu, Mv, Mw, α + 1, β, δ) +
              a[3] * moment_uv_3D1F(Mu, Mv, Mw, α, β + 1, δ) +
              a[4] * moment_uv_3D1F(Mu, Mv, Mw, α, β, δ + 1) +
              0.5 * a[5] * moment_uv_3D1F(Mu, Mv, Mw, α + 2, β, δ) +
              0.5 * a[5] * moment_uv_3D1F(Mu, Mv, Mw, α, β + 2, δ) +
              0.5 * a[5] * moment_uv_3D1F(Mu, Mv, Mw, α, β, δ + 2)
end

