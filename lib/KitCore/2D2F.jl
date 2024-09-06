function discrete_maxwell_2D2F(
    u::AbstractVector{T},
    v::AbstractVector,
    prim::AbstractVector,
    K::Real,
) where {T}
    F = Matrix{T}(undef, length(u), 2)
    @inbounds @. F[:, 1] =
        prim[1] *
        (prim[5] / π) *
        exp(-prim[5] * ((u - prim[2])^2 + (v - prim[3])^2 + (w - prim[4])^2))
    @inbounds @. F[:, 2] = @view(F[:, 1]) * K / (2.0 * prim[4])
    return F
end
function shakhov_part_2D2F(
    u::AbstractVector{T},
    v::AbstractVector,
    H::AbstractVector,
    B::AbstractVector,
    prim::AbstractVector,
    qf::AbstractVector,
    Pr::Real,
    K::Real,
) where {T}
    F⁺ = Matrix{T}(undef, length(u), 2)
    @inbounds @. F⁺[:, 1] =
        0.8 * (1 - Pr) * prim[4]^2 / prim[1] *
        ((u - prim[2]) * qf[1] + (v - prim[3]) * qf[2]) *
        (2 * prim[4] * ((u - prim[2])^2 + (v - prim[3])^2) + K - 5) *
        H
    @inbounds @. F⁺[:, 2] =
        0.8 * (1 - Pr) * prim[4]^2 / prim[1] *
        ((u - prim[2]) * qf[1] + (v - prim[3]) * qf[2]) *
        (2 * prim[4] * ((u - prim[2])^2 + (v - prim[3])^2) + K - 3) *
        B
    return F⁺
end
function heat_flux_2D2F(
    u::AbstractVector{T},
    v::AbstractVector,
    h::AbstractVector,
    b::AbstractVector,
    prim::AbstractVector,
    weight::AbstractVector,
) where{T}
    q = Vector{T}(undef, 2)
    @inbounds q[1] =
        0.5 * (
            sum(@. weight * (u - prim[2]) * ((u - prim[2])^2 + (v - prim[3])^2) * h) +
            sum(@. weight * (u - prim[2]) * b)
        )
    @inbounds q[2] =
        0.5 * (
            sum(@. weight * (v - prim[3]) * ((u - prim[2])^2 + (v - prim[3])^2) * h) +
            sum(@. weight * (v - prim[3]) * b)
        )
    return q
end
function maxwellian_density_2D2F(
    u::AbstractVector,
    v::AbstractVector,
    h::AbstractVector,
    prim::AbstractVector,
    weight::AbstractVector,
    Θ::AbstractVector,
    vn::AbstractVector,
)
    @inbounds SF = sum(@. weight * vn * h * (1.0 - Θ))
    @inbounds SG =
        prim[4] / π *
        sum(@. weight * vn * exp(-prim[4] * ((u - prim[2])^2 + (v - prim[3])^2)) * Θ)
    return -SF / SG
end
function macro_flux_2D2F(
    u::AbstractVector{T},
    v::AbstractVector,
    h::AbstractVector,
    b::AbstractVector,
    weight::AbstractVector,
    vn::AbstractVector,
)where{T}
    flux = Vector{T}(undef,4)
    @inbounds flux[1] = sum(@. weight * vn * h)
    @inbounds flux[2] = sum(@. weight * vn * u * h)
    @inbounds flux[3] = sum(@. weight * vn * v * h)
    @inbounds flux[4] = 0.5 * sum(@. weight * vn * ((u^2 + v^2) * h + b))
    return flux
end
function micro_to_macro_2D2F(u::AbstractVector{T},v::AbstractVector,h::AbstractVector,b::AbstractVector,weight::AbstractVector) where{T}
    w = Vector{T}(undef,4)
    @inbounds w[1] = sum(@. weight * h)
    @inbounds w[2] = sum(@. weight * u * h)
    @inbounds w[3] = sum(@. weight * v * h)
    @inbounds w[4] = 0.5 * sum(@. weight *((u^2 + v^2) * h + b))
    return w
end
function moment_u_2D2F(U::T,V::Real,λ::Real,n::Integer,m::Integer,K::Real) where{T} # n>=1 calculate <u^n> <u^n>_>0 <u^n>_<0 <v^m> <ξ^4>
    Mu_L = Vector{T}(undef, n + 1)
    Mu = Vector{T}(undef, n + 1)
    Mu_R = Vector{T}(undef, n + 1)
    @inbounds begin
        Mu_L[1] = 0.5 * erfc(-√(λ) * U)
        Mu_L[2] =
            U * Mu_L[1] + 0.5 * exp(-λ * U^2) / (√(π * λ))
		Mu_R[1] = 0.5 * erfc(√(λ) * U)
		Mu_R[2] =
		    U * Mu_R[1] - 0.5 * exp(-λ * U^2) / (√(π * λ))
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
		Mξ2 = 0.5 * K / (λ)
		Mξ4 = (K^2 + 2.0 * K) / (4.0 * λ^2)
		Mξ = [one(Mξ2), Mξ2, Mξ4]
	    end
	    return (Mu, Mv, Mξ, Mu_L, Mu_R)
	end

	function moment_uv_2D2F(Mu::AbstractVector{T}, Mv, Mξ, α::Integer, β::Integer, δ::Integer) where{T} # calc <u^αv^βξ^δψ>
	    Muv = Vector{T}(undef, 4)
	    @inbounds begin
		Muv[1] = Mu[α+1] * Mv[β+1] * Mξ[Int(δ / 2 + 1)]
		Muv[2] = Mu[α+2] * Mv[β+1] * Mξ[Int(δ / 2 + 1)]
		Muv[3] = Mu[α+1] * Mv[β+2] * Mξ[Int(δ / 2 + 1)]
		Muv[4] =
		    0.5 * (
			Mu[α+3] * Mv[β+1] * Mξ[Int(δ / 2 + 1)] +
			Mu[α+1] * Mv[β+3] * Mξ[Int(δ / 2 + 1)] +
			Mu[α+1] * Mv[β+1] * Mξ[Int(δ / 2 + 2)]
		    )
	    end
	    return Muv
	end

	function moment_au_2D2F(a, Mu, Mv, Mξ, α::Integer, β::Integer)# calc <au^αv^βψ> 
	    @inbounds a[1] * moment_uv_2D2F(Mu, Mv, Mξ, α, β, 0) +
		      a[2] * moment_uv_2D2F(Mu, Mv, Mξ, α + 1, β, 0) +
		      a[3] * moment_uv_2D2F(Mu, Mv, Mξ, α, β + 1, 0) +
		      0.5 * a[4] * moment_uv_2D2F(Mu, Mv, Mξ, α + 2, β, 0) +
		      0.5 * a[4] * moment_uv_2D2F(Mu, Mv, Mξ, α, β + 2, 0) +
		      0.5 * a[4] * moment_uv_2D2F(Mu, Mv, Mξ, α, β, 2)
	end


