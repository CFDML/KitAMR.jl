abstract type AbstractBoundary end
abstract type AbstractBCType end
abstract type Maxwellian<:AbstractBCType end
const AbstractBC = Union{Vector,Function}
struct Domain{T<:AbstractBCType} <: AbstractBoundary
    id::Integer
    bc::AbstractBC
end
struct Circle{T<:AbstractBCType} <: AbstractBoundary
    center::Vector
    radius::Real
    solid::Bool # If is solid inside the circle, it should be true. Otherwise, it should be false.
    bc::AbstractBC
end

mutable struct SolidCells{DIM,NDF}
    ps_datas::Vector{PS_Data{DIM,NDF}}
    global_midpoints::Vector{Vector} # Global, store all solidcells' midpoints in global field, which is a compromise to the complexity of solidcells' transportation
end
function SolidCells(ps_datas::Vector{PS_Data{DIM,NDF}}) where{DIM,NDF}
    return SolidCells{DIM,NDF}(ps_datas,Vector{Float64}[])
end
mutable struct AuxPoints # Global, one-to-one corresponds to solid cells, intersect point(default)/image point
    midpoint::Vector{Vector}
    mpi_rank::Vector{Integer} # The solidcell corresponding to AuxPoint belongs to which mpi_rank
    # local_index::Vector{Integer} # The solidcell's index on its local processor
end
function AuxPoints()
    return AuxPoints(Vector{Float64}[],Int[])
end
mutable struct GhostIBNodes{DIM,NDF} # Maybe sdf?
    level::AbstractVector
    midpoint::AbstractMatrix
    df::AbstractMatrix
end
const AbstractIBNodes{DIM,NDF} = Union{PS_Data{DIM,NDF},GhostIBNodes{DIM,NDF}}