abstract type AbstractPsData{DIM,NDF} end
abstract type AbstractGhostPsData{DIM,NDF} <: AbstractPsData{DIM,NDF} end
NeighborQuad{DIM,NDF} = Union{AbstractPsData{DIM,NDF}, Nothing}

ChildNum = Union{Val{4},Val{8}}
NeighborNum = Union{Val{2},Val{4}}

BoundaryNeighbor = Val{0}
SameSizeNeighbor = Val{1}
HalfSizeNeighbor = Union{Val{2},Val{4}}
DoubleSizeNeighbor = Val{-1}

mutable struct Neighbor{DIM,NDF}
    data::Vector{Vector{NeighborQuad{DIM,NDF}}}
    state::Vector{Int}
    vg::Vector{Vector{VelocityGradient}}
    # c2r_temp::Vector{Vector{Matrix{Int}}} # Templates indices for velocity grids that are coarser than corresponding grids in neighbors.
    Neighbor(DIM,NDF) = (n = new{DIM,NDF}();
    n.data = Vector{Vector{NeighborQuad{DIM,NDF}}}(undef, 2*DIM);
    n.state = zeros(Int8, 2*DIM);
    n.vg = [VelocityGradient[] for _ in 1:2*DIM];
    # n.c2r_temp = Vector{Vector{Matrix{Int}}}(undef,2*DIM);
    n)
end
mutable struct SolidNeighbor{DIM,NDF,ID} <:AbstractPsData{DIM,NDF}
    bound_enc::Int
    aux_point::Vector{Float64}
    normal::Vector{Float64}
    solid_cell::AbstractPsData{DIM,NDF}
    midpoint::Vector{Float64}
    w::Vector{Float64}
    cvc::CuttedVelocityCells
    vs_data::VS_Data{DIM,NDF}
end
mutable struct PS_Data{DIM,NDF} <: AbstractPsData{DIM,NDF}
    quadid::Cint # The unique identification of the ps_data, is currently used for SolidCells and IB nodes' partition. Only need to be updated before partition. Can be negative for SolidCells/IB nodes for the convenience of boundary_flag
    bound_enc::Int # 0:fluid_cell;>0:bound_enc th boundary's IB_cell;<0: bound_enc th boundary's solid_cell;
    solid_cell_index::Vector{Int}
    ds::Vector{Float64} # DIM
    midpoint::Vector{Float64}
    qf::Vector{Float64}
    w::Vector{Float64}
    sw::Matrix{Float64} # DIM + 2 x DIM
    # ssw::Array{Float64}
    prim::Vector{Float64}
    flux::Vector{Float64}
    vs_data::VS_Data{DIM,NDF}
    neighbor::Neighbor{DIM,NDF}
    PS_Data(DIM,NDF) = (n = new{DIM,NDF}();
    n.quadid = 0;
    n.bound_enc = 0;
    n.solid_cell_index = zeros(SOLID_CELL_ID_NUM);
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
    n.quadid = 0;
    n.bound_enc = 0;
    n.solid_cell_index = zeros(SOLID_CELL_ID_NUM);
    n.ds = ds;
    n.midpoint = midpoint;
    n.w = w;
    n.sw = sw;
    n.flux = zeros(DIM+2);
    n.vs_data = vs_data;
    n)
    PS_Data(DIM,NDF,w, vs_data) = (n = new{DIM,NDF}();
    n.quadid = 0;
    n.bound_enc = 0;
    n.solid_cell_index = zeros(SOLID_CELL_ID_NUM);
    n.w = w;
    n.neighbor = Neighbor(DIM,NDF);
    n.sw = zeros(DIM+2, DIM);
    n.qf = zeros(DIM);
    n.prim = zeros(DIM+2);
    n.flux = zeros(DIM+2);
    n.vs_data = vs_data;
    n)
    PS_Data(DIM,NDF,encs::Vector{Int},w, vs_data) = (n = new{DIM,NDF}();
    n.quadid = 0;
    n.bound_enc = encs[1];
    n.solid_cell_index = encs[2:SOLID_CELL_ID_NUM+1];
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
    bound_enc::Int
    ds::Vector{Cdouble}
    midpoint::Vector{Cdouble}
    w::Vector{Cdouble}
    sw::Matrix{Cdouble}
    vs_data::Ghost_VS_Data{DIM,NDF}
end
# mutable struct SingularSolidNeighbor{DIM,NDF,ID} <: AbstractPsData{DIM,NDF}
#     aux_point::Vector{Float64}
#     normal::Vector{Float64}
#     ghost_ps_data::Ghost_PS_Data{DIM,NDF}
# end
mutable struct GhostInsideSolidData{DIM,NDF} <: AbstractGhostPsData{DIM,NDF} end

abstract type AbstractFace end
abstract type InnerFace <: AbstractFace end
abstract type BoundaryFace <: AbstractFace end
# abstract type AbstractInsideSolidData{DIM,NDF} <: AbstractPsData{DIM,NDF} end
# HangingQuads = Union{Vector{AbstractPsData{DIM,NDF}}, Nothing} where{DIM,NDF}
# struct MissingHangingQuad{DIM,NDF} <: AbstractPsData{DIM,NDF} end
mutable struct InsideSolidData{DIM,NDF} <: AbstractPsData{DIM,NDF} end
AbstractInsideSolidData = Union{InsideSolidData,GhostInsideSolidData}
# struct Face{FT<:AbstractFaceType,T<:HangingQuads}
#     data::PS_Data
#     faceid::Integer
#     hanging_data::T # Is a vector when full is ghost, whose length for 2D = 1 and for 3D = 3.
# end
# function Face(::Type{FT},data::PS_Data,faceid::Integer,hanging_data::T) where{FT<:AbstractFaceType,T<:HangingQuads}
#     return Face{FT,T}(data,faceid,hanging_data)
# end
struct FullFace{DIM,NDF}<:InnerFace
    rot::Float64
    direction::Int
    midpoint::Vector{Float64}
    here_data::PS_Data{DIM,NDF}
    there_data::AbstractPsData{DIM,NDF}
end

struct HangingFace{DIM,NDF}<:InnerFace
    rot::Float64
    direction::Int
    midpoint::Vector{Vector{Float64}}
    here_data::PS_Data{DIM,NDF}
    there_data::Vector{AbstractPsData{DIM,NDF}}
end

struct BackHangingFace{DIM,NDF}<:InnerFace
    rot::Float64
    direction::Int
    midpoint::Vector{Vector{Float64}}
    here_data::Vector{PS_Data{DIM,NDF}}
    there_data::AbstractPsData{DIM,NDF}
end

struct Flux_Data{T<:InnerFace}
    rot::Float64
    direction::Int # face normal direction
    midpoint::Vector{Float64} # midpoint of the face
    here_data::AbstractPsData
    there_data::AbstractPsData
end