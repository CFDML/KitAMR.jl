abstract type AbstractGas end

@with_kw mutable struct Gas <: AbstractGas
    Kn::Float64 = 0.05
    Ma::Float64 = 0.0
    Pr::Float64 = 2/3
    K::Float64 = 1.0 # for 1-D, K = N+2; for 2-D, K = N+1; for 3-D, K = N; where, N denotes the degree of internal freedom
    γ::Float64 = 5/3 # γ=Cₚ/Cᵥ=((N+3)+2)/(N+3). For monotomic gas, γ=5/3.
    ω::Float64 = 0.81
    αᵣ::Float64 = 1.0
    ωᵣ::Float64 = 0.5
    μᵣ::Float64 = ref_vhs_vis(Kn,αᵣ,ωᵣ)
end