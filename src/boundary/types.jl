abstract type AbstractBoundary end
abstract type AbstractBoundaryType end
abstract type Maxwellian<:AbstractBoundaryType end
abstract type SuperSonicInflow <: AbstractBoundaryType end
abstract type SuperSonicOutflow <: AbstractBoundaryType end
abstract type UniformOutflow <: AbstractBoundaryType end
abstract type InterpolatedOutflow <: AbstractBoundaryType end
abstract type AxisSymmetric <: AbstractBoundaryType end
abstract type Period <: AbstractBoundaryType end
const AbstractBCType = Union{Vector,Function}


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
struct Sphere{T<:AbstractBoundaryType} <: AbstractBoundary # 3D circle
    center::Vector
    radius::Real
    solid::Bool
    search_coeffi::Real
    bc::AbstractBCType
    search_radius::Real
    Sphere(::Type{T},center::Vector,radius,solid,search_coeffi,bc) where{T<:AbstractBoundaryType} = new{T}(center,radius,solid,search_coeffi,bc)
    Sphere(c::Sphere{T},ds::Float64) where{T<:AbstractBoundaryType} = new{T}(c.center,
        c.radius,
        c.solid,
        c.search_coeffi,
        c.bc,
        c.search_coeffi*ds
    )
end
const AbstractCircle = Union{Circle,Sphere}

struct Vertices{DIM,T<:AbstractBoundaryType} <: AbstractBoundary
    vertices::Vector{Vector{Float64}} # Vertices of the boundary, sorted in clockwise or counterclockwise order. 
    solid::Bool # Is solid inside the boundary?
    refine_coeffi::Real
    bc::AbstractBCType
    box::Vector{Vector{Float64}} # [[xmin,ymin,zmin],[xmax,ymax,zmax]]
    refine_radius::Real
    Vertices(::Type{T},file::String,solid,refine_coeffi,bc) where{T<:AbstractBoundaryType}= (
        s = CSV.read(file,DataFrame;header=true);
        DIM = length(names(s));
        vertices = DIM==2 ? [[s.x[i],s.y[i]] for i in eachindex(s.x)] : [[s.x[i],s.y[i],s.z[i]] for i in eachindex(s.x)];
        if vertices[1]==vertices[end]
            vertices = vertices[1:end-1]
        end;
        box = DIM==2 ? [[minimum(s.x),minimum(s.y)],[maximum(s.x),maximum(s.y)]] : [[minimum(s.x),minimum(s.y),minimum(s.z)],[maximum(s.x),maximum(s.y),maximum(s.z)]];
        new{DIM,T}(vertices,solid,refine_coeffi,bc,box)
    )
    Vertices(v::Vertices{DIM,T},config) where{DIM,T<:AbstractBoundaryType} =(
        ds_max = maximum([(config[:geometry][2i]-config[:geometry][2i-1])/config[:trees_num][i] for i in 1:config[:DIM]]);
        ds = norm([(config[:geometry][2i]-config[:geometry][2i-1])/config[:trees_num][i]/2^config[:AMR_PS_MAXLEVEL] for i in 1:config[:DIM]]);
        new{DIM,T}(v.vertices,v.solid,v.refine_coeffi,v.bc,[v.box[1].-ds_max,v.box[2].+ds_max],v.refine_coeffi*ds)
    )
end
mutable struct SolidCells{DIM,NDF}
    ps_datas::Vector{PS_Data{DIM,NDF}}
    quadids::Vector{Cint} # A better choice: store all quadid of SolidCells to avoid extensive comparing iteration. Are only needed to be updated before partition.
    #=
    quadids provide the information of solidcells on other processors. Only after this, the exchange of IB nodes can be executeable.
    =#
end
function SolidCells(ps_datas::Vector{PS_Data{DIM,NDF}}) where{DIM,NDF}
    quadids = [x.quadid for x in ps_datas]
    return SolidCells{DIM,NDF}(ps_datas,quadids)
end
mutable struct GhostIBVSData{DIM,NDF} <: AbstractVsData{DIM,NDF}
    vs_num::Int
    level::AbstractVector{Int8}
    weight::AbstractVector{Float64}
    midpoint::AbstractMatrix{Float64}
    df::AbstractMatrix{Float64}
end
mutable struct GhostIBNode{DIM,NDF} # Maybe sdf?
    solid_cell_index::Vector{Int}
    midpoint::Vector{Float64}
    prim::Vector{Float64}
    vs_data::GhostIBVSData{DIM,NDF}
end
const AbstractIBNodes{DIM,NDF} = Union{PS_Data{DIM,NDF},GhostIBNode{DIM,NDF}}
mutable struct IBCells{DIM,NDF}
    IB_nodes::Vector{Vector{AbstractIBNodes{DIM,NDF}}} # SolidCells{IBNodes{}}
    templates::Vector{Vector{Vector{Int}}} # Available templates indices
end
function IBCells(IB_nodes::Vector{Vector{AbstractIBNodes{DIM,NDF}}}) where{DIM,NDF}
    return IBCells{DIM,NDF}(IB_nodes,Vector{Vector{Vector{Int}}}(undef,length(IB_nodes)))
end
mutable struct IBTransferData
    solid_cell_indices::Matrix{Int}
    midpoints::Matrix{Float64}
    prims::Matrix{Float64}
    level::Vector{Int8}
    vs_midpoint::Matrix{Float64}
    df::Matrix{Float64}
end
mutable struct IBBuffer
    sdata::Vector{IBTransferData} # Datas of IB nodes to be sent
    rdata::Vector{IBTransferData} # Datas of IB nodes to be received
    r_vs_nums::Vector{Vector{Int}}
    ghost_nodes::Vector{Vector}
    IBBuffer() = new()
end