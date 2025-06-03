abstract type AbstractVsData{DIM,NDF} end

mutable struct VS_Data{DIM,NDF} <: AbstractVsData{DIM,NDF}
    vs_num::Int
    level::Vector{Int8} # vs_num
    weight::Vector{Float64} # vs_num
    midpoint::Matrix{Float64} # vs_num x DIM
    df::Matrix{Float64} # vs_num x NDF
    sdf::Array{Float64,3} # vs_num x NDF x DIM
    flux::Matrix{Float64} # vs_num x NDF
end

mutable struct Ghost_VS_Data{DIM,NDF} <: AbstractVsData{DIM,NDF}
    vs_num::Int
    level::Vector{Int8} # vs_num
    weight::Vector{Float64} # vs_num
    midpoint::Matrix{Float64} # vs_num x DIM
    df::Matrix{Float64} # vs_num x NDF
    sdf::Array{Float64,3} # vs_num x NDF x DIM
end

# mutable struct Face_VS_Data{DIM,NDF} <: AbstractVsData{DIM,NDF}
#     weight::AbstractVector
#     midpoint::AbstractMatrix
#     vn::AbstractVector
#     df::AbstractMatrix
#     sdf::AbstractMatrix
# end
struct Face_VS_Data{DIM,NDF} # different sides of the face combining the face-velocity-space
    heavi::Vector{Bool}
    weight::AbstractVector{Float64}
    midpoint::AbstractMatrix{Float64}
    vn::AbstractVector{Float64}
    df::AbstractMatrix{Float64}
    sdf::AbstractArray{Float64}
end
function Face_VS_Data(fvd::Face_VS_Data{DIM,NDF},df::AbstractMatrix) where{DIM,NDF}
    return Face_VS_Data{DIM,NDF}(fvd.heavi,fvd.weight,fvd.midpoint,fvd.vn,df,fvd.sdf)
end
struct VelocityTemplates
    indices::Vector{Int} # Indices of the templates
    Ainv::Matrix{Float64} # The inversion of the coefficient marix for bi-linear interpolation
end
mutable struct CuttedVelocityCells
    indices::Vector{Int}
    weight::Vector{Float64}
    gas_dfs::Matrix{Float64}
    solid_dfs::Matrix{Float64}
    gas_midpoints::Matrix{Float64}
    solid_midpoints::Matrix{Float64}
    gas_weights::Vector{Float64}
    solid_weights::Vector{Float64} # Range from 0 to 1, representing the percent of the solid part.
    # templates::Vector{VelocityTemplates}
end
struct VelocityGradient
    vg_template::Matrix{Int} # The complete templates indices for the update of the slope.
    heavi_template::Matrix{Int} # The corresponding templates in the upwind-half-part for the update of the microflux.
end
const AbstractVelocityGradient = Union{VelocityGradient,Nothing}