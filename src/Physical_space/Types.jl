NeighborQuad{DIM,NDF} = Union{AbstractPsData{DIM,NDF},Nothing}

ChildNum = Union{Val{4},Val{8}}
NeighborNum = Union{Val{2},Val{4}}

BoundaryNeighbor = Val{0}
SameSizeNeighbor = Val{1}
HalfSizeNeighbor = Union{Val{2},Val{4}}
DoubleSizeNeighbor = Val{-1}

"""
$(TYPEDEF)
$(TYPEDFIELDS)
"""
mutable struct Neighbor{DIM,NDF}
    data::Vector{Vector{NeighborQuad{DIM,NDF}}}
    state::Vector{Int}
    Neighbor(DIM, NDF) = (
        n = new{DIM,NDF}();
        n.data = Vector{Vector{NeighborQuad{DIM,NDF}}}(undef, 2*DIM);
        n.state = zeros(Int8, 2*DIM);
        n
    )
end

"""
$(TYPEDEF)
$(TYPEDFIELDS)
"""
mutable struct SolidNeighbor{DIM,NDF} <: AbstractPsData{DIM,NDF}
    bound_enc::Int
    """
    Index of the face of the donor cell corresponding to this SolidNeighbor. # The faceid through which the neighbor is the solidneighbor.
    """
    faceid::Int # The faceid through which the neighbor is the solidneighbor.
    """
    Deprecated. # corner neighbor id for corner ghost cells' average.
    """
    av_id::Int # corner neighbor id for corner ghost cells' average.
    """
    Intersect point at the boundary corresponding to this SolidNeighbor.
    """
    aux_point::Vector{Float64}
    """
    Normal vector of the boundary at `aux_point`.
    """
    normal::Vector{Float64}
    """
    Solid_cell corresponding to this SolidNeighbor.
    """
    solid_cell::AbstractPsData{DIM,NDF}
    midpoint::Vector{Float64}
    ds::Vector{Float64}
    w::Vector{Float64}
    sw::Matrix{Float64}
    """
    Struct containing information of cut cell in velocity space.
    """
    cvc::CuttedVelocityCells
    vs_data::VS_Data{DIM,NDF}
end
"""
$(TYPEDEF)
$(TYPEDFIELDS)
"""
mutable struct PS_Data{DIM,NDF} <: AbstractPsData{DIM,NDF}
    """
    Index of the quadrant corresponding to `PS_Data`. # The unique identification of the ps_data, is currently used for SolidCells and IB nodes' partition. Only need to be updated before partition. Can be negative for SolidCells/IB nodes for the convenience of boundary_flag
    """
    quadid::Cint # The unique identification of the ps_data, is currently used for SolidCells and IB nodes' partition. Only need to be updated before partition. Can be negative for SolidCells/IB nodes for the convenience of boundary_flag
    """
    Encoding of the cell type. `0`:fluid cell;`>0`: donor cell of `bound_enc` th immersed boundary;`<0`: solid cell (or solid neighbor) of `bound_enc` th boundary. # 0:fluid_cell;>0:bound_enc th boundary's IB_cell;<0: bound_enc th boundary's solid_cell;
    """
    bound_enc::Int # 0:fluid_cell;>0:bound_enc th boundary's IB_cell;<0: bound_enc th boundary's solid_cell;
    """
    Deprecated.
    """
    solid_cell_index::Vector{Int}
    """
    Faces area (length). # DIM
    """
    ds::Vector{Float64} # DIM
    """
    Coordinates of the cell center.
    """
    midpoint::Vector{Float64}
    """
    Heat flux.
    """
    qf::Vector{Float64}
    """
    Conserved variables.
    """
    w::Vector{Float64}
    """
    Gradients of conserved variables. # DIM + 2 x DIM
    """
    sw::Matrix{Float64} # DIM + 2 x DIM
    """
    Primary variables.
    """
    prim::Vector{Float64}
    """
    Conserved flux.
    """
    flux::Vector{Float64}
    """
    Data of the velocity space.
    """
    vs_data::VS_Data{DIM,NDF}
    """
    Neighbor of the cell.
    """
    neighbor::Neighbor{DIM,NDF}
    PS_Data(DIM, NDF; kwargs...) = (
        n = new{DIM,NDF}();
        n.quadid = haskey(kwargs, :quadid) ? kwargs[:quadid] : 0;
        n.bound_enc = haskey(kwargs, :bound_enc) ? kwargs[:bound_enc] : 0;
        n.solid_cell_index = haskey(kwargs, :solid_cell_index) ?
                             kwargs[:solid_cell_index] : zeros(SOLID_CELL_ID_NUM);
        n.ds = haskey(kwargs, :ds) ? kwargs[:ds] : zeros(DIM);
        n.midpoint = haskey(kwargs, :midpoint) ? kwargs[:midpoint] : zeros(DIM);
        n.qf = haskey(kwargs, :qf) ? kwargs[:qf] : zeros(DIM);
        n.w = haskey(kwargs, :w) ? kwargs[:w] : zeros(DIM+2);
        n.sw = haskey(kwargs, :sw) ? kwargs[:sw] : zeros(DIM+2, DIM);
        n.prim = haskey(kwargs, :prim) ? kwargs[:prim] : zeros(DIM+2);
        n.flux = haskey(kwargs, :flux) ? kwargs[:flux] : zeros(DIM+2);
        if haskey(kwargs, :vs_data)
            n.vs_data = kwargs[:vs_data]
        end;
        n.neighbor = haskey(kwargs, :neighbor) ? kwargs[:neighbor] : Neighbor(DIM, NDF);
        n
    )
    PS_Data(ps_data::PS_Data{DIM,NDF}; kwargs...) where {DIM,NDF} = (
        n = new{DIM,NDF}();
        for field in fieldnames(PS_Data{DIM,NDF})
            setfield!(
                n,
                field,
                haskey(kwargs, field) ? kwargs[field] : getfield(ps_data, field),
            )
        end;
        n
    )
end

"""
$(TYPEDEF)
$(TYPEDFIELDS)
"""
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

"""
$(TYPEDEF)
"""
mutable struct GhostInsideSolidData{DIM,NDF} <: AbstractGhostPsData{DIM,NDF} end

"""
$(TYPEDEF)
"""
mutable struct InsideSolidData{DIM,NDF} <: AbstractPsData{DIM,NDF} end
AbstractInsideSolidData = Union{InsideSolidData,GhostInsideSolidData}

"""
$(TYPEDEF)
Face shared by cells with the same refinement level.
$(TYPEDFIELDS)
"""
struct FullFace{DIM,NDF}<:InnerFace
    """
    Rotational coefficient of the face. `-1` if `here_data` locates at the negative direction relative to the face.
    """
    rot::Float64
    """
    Direction of the face's normal. 1, 2, 3 correspond to x, y, z respectively.
    """
    direction::Int
    """
    Coordinates of the face's center.
    """
    midpoint::Vector{Float64}
    """
    `PS_Data` on one side of the face. It is guaranteed to be local.
    """
    here_data::PS_Data{DIM,NDF}
    """
    `AbstractPsData` on the other side of the face. It can be a ghost one.
    """
    there_data::AbstractPsData{DIM,NDF}
end

"""
$(TYPEDEF)
Faces shared by cells with different refinement levels. The cell with lower level is local.
$(TYPEDFIELDS)
"""
struct HangingFace{DIM,NDF}<:InnerFace
    rot::Float64
    direction::Int
    midpoint::Vector{Vector{Float64}}
    here_data::PS_Data{DIM,NDF}
    """
    `AbstractPsData`s with higher refinement level.
    """
    there_data::Vector{AbstractPsData{DIM,NDF}}
end

"""
$(TYPEDEF)
Faces shared by cells with different refinement levels. The cell with lower level is ghost.
$(TYPEDFIELDS)
"""
struct BackHangingFace{DIM,NDF}<:InnerFace
    rot::Float64
    direction::Int
    midpoint::Vector{Vector{Float64}}
    here_data::Vector{PS_Data{DIM,NDF}}
    there_data::AbstractPsData{DIM,NDF}
end

"""
$(TYPEDEF)
Struct for flux calculation corresponding to a single face.
$(TYPEDFIELDS)
"""
struct FluxData{T<:InnerFace}
    rot::Float64
    direction::Int # face normal direction
    midpoint::Vector{Float64} # midpoint of the face
    here_data::AbstractPsData
    there_data::AbstractPsData
end

mutable struct MeshData
    in_box::Int
    in_solid::Bool
    in_search_radius::Int # In i-th ib's search radius. If not, i is set to 0.
    is_ghost_cell::Bool
end
function MeshData()
    return MeshData(0, false, 0, false)
end
