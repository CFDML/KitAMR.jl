"""
$(TYPEDEF)
$(TYPEDFIELDS)
"""
mutable struct Gas <: AbstractGas
    "Default is `0.05`."
    Kn::Float64
    "Default is `2/3`."
    Pr::Float64
    "The translational internal DOF.Default is `1.0`."
    K::Float64 # for 1-D, K = N+2; for 2-D, K = N+1; for 3-D, K = N; where, N denotes the degree of internal freedom (containing rotational and viberational, etc., but not translational.)
    "Default is `5/3`."
    γ::Float64 # γ=Cₚ/Cᵥ=((N+3)+2)/(N+3). For monotomic gas, N=0 and γ=5/3.
    "Default is `0.5`."
    ω::Float64
    "Default is `1.0`."
    αᵣ::Float64
    "Default is `0.81`."
    ωᵣ::Float64
    "The reference viscosity. Default is computed by [`ref_vhs_vis`](@ref)."
    μᵣ::Float64
    "The reference temperature of the molecular model. Default is `1.0`."
    T_ref::Float64 # The reference temperature of the molecular model
end
function Gas(;kwargs...)
    Kn = haskey(kwargs,:Kn) ? kwargs[:Kn] : 0.05
    Pr = haskey(kwargs,:Pr) ? kwargs[:Pr] : 2/3
    K = haskey(kwargs,:K) ? kwargs[:K] : 1.0
    γ = haskey(kwargs,:γ) ? kwargs[:γ] : 5/3
    ω = haskey(kwargs,:ω) ? kwargs[:ω] : 0.5
    αᵣ = haskey(kwargs,:αᵣ) ? kwargs[:αᵣ] : 1.0
    ωᵣ = haskey(kwargs,:ωᵣ) ? kwargs[:ωᵣ] : 0.81
    μᵣ = haskey(kwargs,:μᵣ) ? kwargs[:μᵣ] : ref_vhs_vis(Kn,αᵣ,ωᵣ)
    T_ref = haskey(kwargs,:T_ref) ? kwargs[:T_ref] : 1.0
    return Gas(
        Kn,Pr,K,γ,ω,αᵣ,ωᵣ,μᵣ,T_ref
    )
end