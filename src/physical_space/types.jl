abstract type AbstractPsData{DIM,NDF} end
abstract type AbstractGhostPsData{DIM,NDF} <: AbstractPsData{DIM,NDF} end
Neighbor_Quad{DIM,NDF} = Union{AbstractPsData{DIM,NDF}, Nothing}

ChildNum = Union{Val{4},Val{8}}
NeighborNum = Union{Val{2},Val{4}}

BoundaryNeighbor = Val{0}
SameSizeNeighbor = Val{1}
HalfSizeNeighbor = Union{Val{2},Val{4}}
DoubleSizeNeighbor = Val{-1}

mutable struct Neighbor{DIM,NDF}
    data::Vector{Vector{Neighbor_Quad{DIM,NDF}}}
    state::Vector{Int}
    Neighbor(DIM,NDF) = (n = new{DIM,NDF}();
    n.data = Vector{Vector{Neighbor_Quad{DIM,NDF}}}(undef, 2*DIM);
    n.state = zeros(Int8, 2*DIM);
    n)
end

mutable struct PS_Data{DIM,NDF} <: AbstractPsData{DIM,NDF}
    quadid::Cint # The unique identification of the ps_data, is currently used for SolidCells and IB nodes' partition. Only need to be updated before partition. Can be negative for SolidCells/IB nodes for the convenience of boundary_flag
    ds::Vector{Float64} # DIM
    midpoint::Vector{Float64}
    qf::Vector{Float64}
    w::Vector{Float64}
    sw::Matrix{Float64} # DIM + 2 x DIM
    prim::Vector{Float64}
    flux::Vector{Float64}
    vs_data::VS_Data{DIM,NDF}
    neighbor::Neighbor{DIM,NDF}
    PS_Data(DIM,NDF) = (n = new{DIM,NDF}();
    n.ds = zeros(DIM);
    n.midpoint = zeros(DIM);
    n.qf = zeros(DIM);
    n.w = zeros(DIM+2);
    n.sw = zeros(DIM+2, DIM);
    n.prim = zeros(DIM+2);
    n.flux = zeros(DIM+2);
    n.neighbor = Neighbor(DIM,NDF);
    n)
    PS_Data(DIM,NDF,ds, midpoint, w, sw, vs_data) = (n = new{DIM,NDF}();
    n.ds = ds;
    n.midpoint = midpoint;
    n.w = w;
    n.sw = sw;
    n.flux = zeros(DIM+2);
    n.vs_data = vs_data;
    n)
    PS_Data(DIM,NDF,w, vs_data) = (n = new{DIM,NDF}();
    n.w = w;
    n.neighbor = Neighbor(DIM,NDF);
    n.sw = zeros(DIM+2, DIM);
    n.qf = zeros(DIM);
    n.prim = zeros(DIM+2);
    n.flux = zeros(DIM+2);
    n.vs_data = vs_data;
    n)
end

mutable struct Ghost_PS_Data{DIM,NDF}<:AbstractGhostPsData{DIM,NDF}
    ds::Vector{Cdouble}
    midpoint::Vector{Cdouble}
    w::Vector{Cdouble}
    sw::Matrix{Cdouble}
    vs_data::Ghost_VS_Data{DIM,NDF}
end
mutable struct GhostInsideSolidData{DIM,NDF} <: AbstractGhostPsData{DIM,NDF} end

abstract type AbstractFaceType end
abstract type InnerFace <: AbstractFaceType end
abstract type BoundaryFace <: AbstractFaceType end
# abstract type AbstractInsideSolidData{DIM,NDF} <: AbstractPsData{DIM,NDF} end
HangingQuads = Union{Vector{AbstractPsData{DIM,NDF}}, Nothing} where{DIM,NDF}
struct MissingHangingQuad{DIM,NDF} <: AbstractPsData{DIM,NDF} end
mutable struct InsideSolidData{DIM,NDF} <: AbstractPsData{DIM,NDF} end
AbstractInsideSolidData = Union{InsideSolidData,GhostInsideSolidData}
struct Face{FT<:AbstractFaceType,T<:HangingQuads}
    data::PS_Data
    faceid::Integer
    hanging_data::T # Is a vector when full is ghost, whose length for 2D = 1 and for 3D = 3.
end
function Face(::Type{FT},data::PS_Data,faceid::Integer,hanging_data::T) where{FT<:AbstractFaceType,T<:HangingQuads}
    return Face{FT,T}(data,faceid,hanging_data)
end
