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
    search_radius::Real # Search radius determines the max distance between IB nodes and aux_points. Default is 
    bc::AbstractBC
end

mutable struct SolidCells{DIM,NDF}
    ps_datas::Vector{PS_Data{DIM,NDF}}
    # global_midpoints::Vector{Vector} # Global, store all solidcells' midpoints in global field, which is a compromise to the complexity of solidcells' transportation
    qudids::Vector{Cint} # A better choice: store all quadid of SolidCells to avoid extensive comparing iteration. Are only needed to be updated before partition.
end
function SolidCells(ps_datas::Vector{PS_Data{DIM,NDF}}) where{DIM,NDF}
    quadids = [x.quadid for x in ps_datas]
    return SolidCells{DIM,NDF}(ps_datas,quadids)
end
# mutable struct AuxPoints # Global, one-to-one corresponds to solid cells, intersect point(default)/image point
#     midpoint::Vector{Vector}
#     mpi_rank::Vector{Integer} # The solidcell corresponding to AuxPoint belongs to which mpi_rank
#     # local_index::Vector{Integer} # The solidcell's index on its local processor
# end

mutable struct GhostIBNodes{DIM,NDF} # Maybe sdf?
    midpoint::Vector{Float64}
    level::AbstractVector{Int8}
    df::AbstractMatrix
end
const AbstractIBNodes{DIM,NDF} = Union{PS_Data{DIM,NDF},GhostIBNodes{DIM,NDF}}
mutable struct IBCells
    IB_nodes::Vector{Vector{AbstractIBNodes}} # boundaries{SolidCells{IBNodes{}}}
    quadids::Vector{Vector{Cint}} # Global quadids. Are only needed to be updated before partition.
end