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
    Neighbor(DIM,NDF) = (n = new{DIM,NDF}();
    n.data = Vector{Vector{NeighborQuad{DIM,NDF}}}(undef, 2*DIM);
    n.state = zeros(Int8, 2*DIM);
    n)
end
mutable struct SolidNeighbor{DIM,NDF} <:AbstractPsData{DIM,NDF}
    bound_enc::Int
    faceid::Int # The faceid through which the neighbor is the solidneighbor.
    av_id::Int # corner neighbor id for corner ghost cells' average.
    aux_point::Vector{Float64}
    normal::Vector{Float64}
    solid_cell::AbstractPsData{DIM,NDF}
    midpoint::Vector{Float64}
    ds::Vector{Float64}
    w::Vector{Float64}
    sw::Matrix{Float64}
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
    prim::Vector{Float64}
    flux::Vector{Float64}
    vs_data::VS_Data{DIM,NDF}
    neighbor::Neighbor{DIM,NDF}
    PS_Data(DIM,NDF;kwargs...)=(n = new{DIM,NDF}();
        n.quadid = haskey(kwargs,:quadid) ? kwargs[:quadid] : 0;
        n.bound_enc = haskey(kwargs,:bound_enc) ? kwargs[:bound_enc] : 0;
        n.solid_cell_index = haskey(kwargs,:solid_cell_index) ? kwargs[:solid_cell_index] : zeros(SOLID_CELL_ID_NUM);
        n.ds = haskey(kwargs,:ds) ? kwargs[:ds] : zeros(DIM);
        n.midpoint = haskey(kwargs,:midpoint) ? kwargs[:midpoint] : zeros(DIM);
        n.qf = haskey(kwargs,:qf) ? kwargs[:qf] : zeros(DIM);
        n.w = haskey(kwargs,:w) ? kwargs[:w] : zeros(DIM+2);
        n.sw = haskey(kwargs,:sw) ? kwargs[:sw] : zeros(DIM+2, DIM);
        n.prim = haskey(kwargs,:prim) ? kwargs[:prim] : zeros(DIM+2);
        n.flux = haskey(kwargs,:flux) ? kwargs[:flux] : zeros(DIM+2);
        if haskey(kwargs,:vs_data)
           n.vs_data = kwargs[:vs_data]
        end;
        n.neighbor = haskey(kwargs,:neighbor) ? kwargs[:neighbor] : Neighbor(DIM,NDF);
        n
    )
    PS_Data(ps_data::PS_Data{DIM,NDF};kwargs...) where{DIM,NDF}=(n = new{DIM,NDF}();
        for field in fieldnames(PS_Data{DIM,NDF})
            setfield!(n,field,haskey(kwargs,field) ? kwargs[field] : getfield(ps_data,field))
        end;
        n
    )
end
mutable struct Ghost_PS_Data{DIM,NDF}<:AbstractGhostPsData{DIM,NDF}
    owner_rank::Int
    quadid::Int # Not continuous! Only provide an order among ghost_datas from the same owner rank.
    bound_enc::Int
    ds::Vector{Cdouble}
    midpoint::Vector{Cdouble}
    w::Vector{Cdouble}
    sw::Matrix{Cdouble}
    vs_data::Ghost_VS_Data{DIM,NDF}
end
mutable struct GhostInsideSolidData{DIM,NDF} <: AbstractGhostPsData{DIM,NDF} end


mutable struct InsideSolidData{DIM,NDF} <: AbstractPsData{DIM,NDF} end
AbstractInsideSolidData = Union{InsideSolidData,GhostInsideSolidData}
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