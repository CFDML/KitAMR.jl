abstract type AbstractVsData{DIM,NDF} end

mutable struct VS_Data{DIM,NDF} <: AbstractVsData{DIM,NDF}
    vs_num::Int
    level::Vector{Int} # vs_num
    weight::Vector{Float64} # vs_num
    midpoint::Matrix{Float64} # vs_num x DIM
    df::Matrix{Float64} # vs_num x NDF
    sdf::Array{Float64,3} # vs_num x NDF x DIM
    flux::Matrix{Float64} # vs_num x NDF
end

mutable struct Ghost_VS_Data{DIM,NDF} <: AbstractVsData{DIM,NDF}
    vs_num::Int
    level::Vector{Int} # vs_num
    weight::Vector{Float64} # vs_num
    midpoint::Matrix{Float64} # vs_num x DIM
    df::Matrix{Float64} # vs_num x NDF
    sdf::Array{Float64,3} # vs_num x NDF x DIM
end

mutable struct Face_VS_Data{DIM,NDF} <: AbstractVsData{DIM,NDF}
    weight::AbstractVector
    midpoint::AbstractMatrix
    vn::AbstractVector
    df::AbstractMatrix
    sdf::AbstractMatrix
end