module KitCore
    using SpecialFunctions
    using ExportAll
    include("2D.jl")
    include("2D2F.jl")
    include("3D.jl")
    include("3D1F.jl")
    function time_int(τ::T, Δt)where{T}
        Mt = Vector{T}(undef, 5)
        @inbounds begin
            Mt[4] = τ * (1.0 - exp(-Δt / τ))
            Mt[5] = -τ * Δt * exp(-Δt / τ) + τ * Mt[4]
            Mt[1] = Δt - Mt[4]
            Mt[2] = -τ * Mt[1] + Mt[5]
            Mt[3] = 0.5 * Δt^2 - τ * Mt[1]
        end
        return Mt
    end
    
    @exportAll()
end # module KitCore
