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

mutable struct GhostIBNodes{DIM,NDF} # Maybe sdf?
    level::AbstractVector
    midpoint::AbstractMatrix
    df::AbstractMatrix
end
const AbstractIBNodes{DIM,NDF} = Union{PS_Data{DIM,NDF},GhostIBNodes{DIM,NDF}}