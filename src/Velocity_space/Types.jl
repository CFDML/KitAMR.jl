"""
$(TYPEDEF)
$(TYPEDFIELDS)
"""
mutable struct VsData{DIM,NDF} <: AbstractVsData{DIM,NDF}
    "Number of velocity cells."
    vs_num::Int
    "Refinement level of velocity cells."
    level::Vector{Int8} # vs_num
    "Quadrature weight of velocity cells."
    weight::Vector{Float64} # vs_num
    "Coordinates of centers of velocity cells in velocity space."
    midpoint::Matrix{Float64} # vs_num x DIM
    "Discretized distribution function."
    df::Matrix{Float64} # vs_num x NDF
    "Spatial gradients (df/dx) of distribution function."
    sdf::Array{Float64,3} # vs_num x NDF x DIM
    "Numerical flux of distribution function."
    flux::Matrix{Float64} # vs_num x NDF
end

"""
$(TYPEDEF)
$(TYPEDFIELDS)
"""
mutable struct Ghost_VsData{DIM,NDF} <: AbstractVsData{DIM,NDF}
    vs_num::Int
    level::Vector{Int8} # vs_num
    weight::Vector{Float64} # vs_num
    midpoint::Matrix{Float64} # vs_num x DIM
    df::Matrix{Float64} # vs_num x NDF
    sdf::Array{Float64,3} # vs_num x NDF x DIM
end

"""
$(TYPEDEF)
Information of the velocity space corresponding to the face during the calculation of numerical flux.
$(TYPEDFIELDS)
"""
struct FaceVsData{DIM,NDF} # different sides of the face combining the face-velocity-space
    "Upwinding flags. Upwinding velocities are set to `true`."
    heavi::Vector{Bool}
    weight::AbstractVector{Float64}
    midpoint::AbstractMatrix{Float64}
    "Normal component of velocities across the face."
    vn::AbstractVector{Float64}
    df::AbstractMatrix{Float64}
    sdf::AbstractArray{Float64}
end
function FaceVsData(fvd::FaceVsData{DIM,NDF},df::AbstractMatrix) where{DIM,NDF}
    return FaceVsData{DIM,NDF}(fvd.heavi,fvd.weight,fvd.midpoint,fvd.vn,df,fvd.sdf)
end
struct VelocityTemplates
    indices::Vector{Int} # Indices of the templates
    Ainv::Matrix{Float64} # The inversion of the coefficient marix for bi-linear interpolation
end

"""
$(TYPEDEF)
Structure for cut cell in velocity space used in immersed boundaries.
$(TYPEDFIELDS)
"""
mutable struct CuttedVelocityCells
    "Indices of cut cells in [`SolidNeighbor`](@ref)'s [`VsData`](@ref)."
    indices::Vector{Int}
    "Quadrature weight of velocity cells. Specially, the weight of cut cells is set to 0."
    weight::Vector{Float64}
    "Distribution function of gas part corresponding to velocity at the center of cut cells."
    gas_dfs::Matrix{Float64}
    "Distribution function of solid part corresponding to velocity at the center of cut cells."
    solid_dfs::Matrix{Float64}
    "Areas of gas part in cut cells."
    gas_weights::Vector{Float64}
    "Areas of solid part in cut cells."
    solid_weights::Vector{Float64}
end
mutable struct VS_Projection{DIM}
    c2r_index::Matrix{Int} #
    c2r_offset::Matrix{Int}
    r2c_index::Matrix{Int}
    r2c_offset::Matrix{Int}
end

struct Gauss_Hermite{NUM} <: AbstractQuadrature
    vcoords::AbstractVector
    weights::AbstractVector
end
function Gauss_Hermite(;NP=28)
    if NP==28
        vcoords = [ -0.5392407922630E+01, -0.4628038787602E+01, -0.3997895360339E+01, -0.3438309154336E+01,
                    -0.2926155234545E+01, -0.2450765117455E+01, -0.2007226518418E+01, -0.1594180474269E+01,
                    -0.1213086106429E+01, -0.8681075880846E+00, -0.5662379126244E+00, -0.3172834649517E+00,
                    -0.1331473976273E+00, -0.2574593750171E-01, +0.2574593750171E-01, +0.1331473976273E+00,
                    +0.3172834649517E+00, +0.5662379126244E+00, +0.8681075880846E+00, +0.1213086106429E+01,
                    +0.1594180474269E+01, +0.2007226518418E+01, +0.2450765117455E+01, +0.2926155234545E+01,
                    +0.3438309154336E+01, +0.3997895360339E+01, +0.4628038787602E+01, +0.5392407922630E+01 ]
        weights = [ +0.2070921821819E-12, +0.3391774320172E-09, +0.6744233894962E-07, +0.3916031412192E-05,
                    +0.9416408715712E-04, +0.1130613659204E-02, +0.7620883072174E-02, +0.3130804321888E-01,
                    +0.8355201801999E-01, +0.1528864568113E+00, +0.2012086859914E+00, +0.1976903952423E+00,
                    +0.1450007948865E+00, +0.6573088665062E-01, +0.6573088665062E-01, +0.1450007948865E+00,
                    +0.1976903952423E+00, +0.2012086859914E+00, +0.1528864568113E+00, +0.8355201801999E-01,
                    +0.3130804321888E-01, +0.7620883072174E-02, +0.1130613659204E-02, +0.9416408715712E-04,
                    +0.3916031412192E-05, +0.6744233894962E-07, +0.3391774320172E-09, +0.2070921821819E-12 ]
        return Gauss_Hermite{NP}(vcoords,weights)
    else
        throw(`Gauss_Hermite not defined yet.`)
    end
end

struct Velocity_Resolution
    density::Float64
    energy::Float64
end

struct Velocity_Adaptive_Data
    vr::Velocity_Resolution
    va_flags::Vector{Bool}
end