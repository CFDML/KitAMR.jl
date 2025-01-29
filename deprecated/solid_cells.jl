"""
The flag represents wheather a cell at `midpoint` in velocity space is closer than `Δ`/2 to the discontinuity caused by the immersed boundary method.
"""
function discontinuity_flag(n::Vector{Float64},midpoint::AbstractVector;Δ) # Δ represents the width of the discontinuity
    dot(midpoint,n)<=Δ/2&&return true
    return false
end
function discontinuity_flag(n::Vector{Float64},midpoints::AbstractMatrix;Δ)
    for midpoint in eachrow(midpoints)
        dot(midpoint,n)<=Δ/2&&return true
    end
    return false
end