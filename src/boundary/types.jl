abstract type AbstractBoundary end
abstract type AbstractBoundaryType end
abstract type Maxwellian<:AbstractBoundaryType end
abstract type SuperSonicInflow <: AbstractBoundaryType end
abstract type SuperSonicOutflow <: AbstractBoundaryType end
abstract type UniformOutflow <: AbstractBoundaryType end
abstract type AxisSymmetric <: AbstractBoundaryType end
abstract type Period <: AbstractBoundaryType end
const AbstractBCType = Union{Vector,Function}
abstract type AbstractIBSortType end
abstract type AbstractIBInterpolateType end

struct DistanceIBSort<:AbstractIBSortType end
struct BilinearIBInterpolate<:AbstractIBInterpolateType end
struct Domain{T<:AbstractBoundaryType} <: AbstractBoundary
    id::Int
    bc::AbstractBCType
    Domain(T,id) = new{T}(id)
    Domain(T,id,bc) = new{T}(id,bc)
end
struct DomainFace{DIM,NDF,T}<:BoundaryFace
    rot::Float64
    direction::Int
    midpoint::Vector{Float64}
    domain::Domain{T}
    ps_data::PS_Data{DIM,NDF}
end
struct Circle{T<:AbstractBoundaryType} <: AbstractBoundary
    center::Vector
    radius::Real
    solid::Bool # If is solid inside the circle, it should be true. Otherwise, it should be false.
    search_coeffi::Real # The ratio of the search radius and the minimal mesh scale, which determines the maximal distance between IB nodes and aux_points.
    bc::AbstractBCType
    search_radius::Real
    Circle(::Type{T},center::Vector,radius,solid,search_coeffi,bc) where{T<:AbstractBoundaryType} = new{T}(center,radius,solid,search_coeffi,bc)
    Circle(c::Circle{T},ds::Float64) where{T<:AbstractBoundaryType} = new{T}(c.center,
        c.radius,
        c.solid,
        c.search_coeffi,
        c.bc,
        c.search_coeffi*ds
    )
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
    weight::AbstractVector{Float64}
    midpoint::AbstractMatrix{Float64}
    df::AbstractMatrix{Float64}
    # sdf::AbstractArray{Float64}
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
    level::Vector{Int8}
    weight::Vector{Float64}
    vs_midpoint::Matrix{Float64}
    df::Matrix{Float64}
    # sdf::Array{Float64}
end
mutable struct IBBuffer
    sdata::Vector{IBTransferData} # Datas of IB nodes to be sent
    rdata::Vector{IBTransferData} # Datas of IB nodes to be received
    r_vs_nums::Vector{Vector{Int}}
    IBBuffer() = new()
end