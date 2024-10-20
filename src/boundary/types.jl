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
    quadids::Vector{Cint} # A better choice: store all quadid of SolidCells to avoid extensive comparing iteration. Are only needed to be updated before partition.
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
mutable struct GhostIBVSData{DIM,NDF} <: AbstractVsData{DIM,NDF}
    vs_num::Int
    level::AbstractVector{Int8}
    df::AbstractMatrix{Float64}
end
mutable struct GhostIBNode{DIM,NDF} # Maybe sdf?
    midpoint::Vector{Float64}
    prim::Vector{Float64}
    vs_data::GhostIBVSData{DIM,NDF}
end
const AbstractIBNodes{DIM,NDF} = Union{PS_Data{DIM,NDF},GhostIBNode{DIM,NDF}}
mutable struct IBCells
    IB_nodes::Vector{Vector{AbstractIBNodes}} # boundaries{SolidCells{IBNodes{}}}
    # quadids::Vector{Vector{Cint}} # Global quadids. Are only needed to be updated before partition.
end
mutable struct IBTransferData
    midpoints::Matrix{Float64}
    prims::Matrix{Float64}
    level::Vector{Int}
    df::Matrix{Float64}
end
mutable struct IBBuffer
    sdata::Vector{IBTransferData} # Datas of IB nodes to be sent
    rdata::Vector{IBTransferData} # Datas of IB nodes to be received
    r_vs_nums::Vector{Vector{Int}}
    IBBuffer() = new()
end