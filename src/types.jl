abstract type AbstractPsData end
abstract type AbstractVsData end
abstract type AbstractFaceWrap end
struct Face{T<:Union{AbstractPsData,Nothing}}
    data::AbstractPsData
    faceid::Int
    bound::Int # 0: inner 1: boundary 2: hanging mirror
    hanging_data::T # Todo: for 3D, hanging_data should be a vector
end
# mutable struct Ghost_Data{T}
#     ds::SVector{DIM,Cdouble}
#     midpoint::SVector{DIM,Cdouble}
#     w::SVector{DIM + 2,Cdouble}
#     local_vs_num::Int
#     micro::T # ds, midpoint, df, sdf
# end
# mutable struct Ghost_Data_Temp
#     ds::Vector{Float64}
#     midpoint::Vector{Float64}
#     w::Vector{Float64}
#     local_vs_num::Int
#     micro::Vector{Float64}
# end
# mutable struct Ghost_Slope{T}
#     sdf::T
# end
# mutable struct Ghost_Slope_Temp
#     sdf::Vector{Float64}
# end
# function Ghost_Slope_T(vs_num::Integer)
#     return Ghost_Slope{SVector{vs_num * NDF * DIM,Cdouble}}
# end
# mutable struct Ghost_VS_Structure{T1,T2,T3}
#     weight::T1
#     level::T2
#     midpoint::T3
# end
# mutable struct Ghost_VS_Structure_Temp
#     weight::Vector{Float64}
#     level::Vector{Int}
#     midpoint::Vector{Float64}
# end
# function Ghost_VS_Structure_T(vs_num::Integer)
#     return Ghost_VS_Structure{SVector{vs_num,Float64},SVector{vs_num,Int},SVector{vs_num*DIM,Float64}}
# end
# mutable struct Ghost_VS_Structure{T1,T2}
#     neighbor_num::NTuple{2*DIM,Int}
#     neighbor_states::T1 # [2*DIM*2^(DIM-1)]*[vs_num]
#     level::T2
#     index::T2
# end
# function Ghost_VS_Structure_T(vs_num::Integer)
#     return Ghost_VS_Structure{NTuple{2 * DIM * vs_num*2^(DIM-1),Int},NTuple{vs_num,Int}}
# end
mutable struct Ghost_VS_Data <: AbstractVsData
    vs_num::Int
    level::Vector{Int} # vs_num
    weight::Vector{Float64} # vs_num
    midpoint::Matrix{Float64} # vs_num x DIM
    df::Matrix{Float64} # vs_num x NDF
    sdf::Array{Float64,3} # vs_num x NDF x DIM
end
# function Ghost_Data_T(vs_num::Integer)
#     return Ghost_Data{SVector{NDF * vs_num,Cdouble}}
# end
mutable struct P4est_PS_Data
    ps_data::Ptr{Nothing}
end
function P4est_PS_Data()
    return P4est_PS_Data(C_NULL)
end
# abstract type AbstractSpaceStructure end
mutable struct Neighbor{T}
    data::Array{Array{T,1},1}
    state::Vector{Int64}
end
# mutable struct VS_Structure <: AbstractSpaceStructure 
#     neighbor_states::Vector{Vector{Vector{Int}}} #  dir{neighbor_num{state{}}};-x, +x, -y, +y, -z, +z; 1, 2^DIM, -1
#     level::Vector{Int} # begin from 1
#     index::Vector{Int} # 1..2^DIM
#     # inherent::Int
# end
mutable struct Ghost_PS_Data{T<:AbstractVsData} <: AbstractPsData
    ds::Vector{Cdouble}
    midpoint::Vector{Cdouble}
    w::Vector{Cdouble}
    vs_data::T
end

#test
function Neighbor(::Type{T}) where {T}
    return Neighbor(Vector{Vector{T}}(undef, 2 * DIM), zeros(Int64, 2 * DIM))
end

mutable struct VS_Data <: AbstractVsData
    vs_num::Int
    level::Vector{Int} # vs_num
    weight::Vector{Float64} # vs_num
    midpoint::Matrix{Float64} # vs_num x DIM
    df::Matrix{Float64} # vs_num x NDF
    sdf::Array{Float64,3} # vs_num x NDF x DIM
    flux::Matrix{Float64} # vs_num x NDF
end

mutable struct PS_Data <: AbstractPsData
    ds::Vector{Cdouble} # DIM
    midpoint::Vector{Float64}
    qf::Vector{Float64}
    w::Vector{Float64}
    sw::Matrix{Float64} # DIM + 2 x DIM
    prim::Vector{Float64}
    flux::Vector{Float64}
    vs_data::VS_Data
    # vs_structure::VS_Structure
    neighbor::Neighbor{Union{AbstractPsData,Nothing}}
    PS_Data() = (n = new();
    n.ds = zeros(DIM);
    n.midpoint = zeros(DIM);
    n.qf = zeros(DIM);
    n.w = zeros(DIM + 2);
    n.sw = zeros(DIM + 2, DIM);
    n.prim = zeros(DIM + 2);
    n.flux = zeros(DIM + 2);
    n.neighbor = Neighbor(Union{AbstractPsData,Nothing});
    n)
    PS_Data(ds, midpoint, w, sw, vs_data) = (n = new();
    n.ds = ds;
    n.midpoint = midpoint;
    n.w = w;
    n.sw = sw;
    n.flux = zeros(DIM + 2);
    n.vs_data = vs_data;
    n)
    PS_Data(w, vs_data) = (n = new();
    n.w = w;
    n.neighbor = Neighbor(Union{AbstractPsData,Nothing});
    n.sw = zeros(DIM + 2, DIM);
    n.qf = zeros(DIM);
    n.prim = zeros(DIM + 2);
    n.flux = zeros(DIM + 2);
    n.vs_data = vs_data;
    n)
end
mutable struct Global_Data{T1,T2,T3,T4,T5,T6,T7}
    geometry::T1
    trees_num::T2
    gas::T3
    ic::T4
    bc::T5
    quadrature::T6
    vs_trees_num::T7
    ghost::Ptr{p4est_ghost_t}
    mesh::Ptr{p4est_mesh_t}
    vs_num::Int
    gradmax::Float64
    adapt_step::Int
    vs_adapt_step::Int
    partition_step::Int
    ps_save_index::Int
    vs_save_index::Int
end
@with_kw mutable struct Gas{T1<:Real} <: AbstractGas
    Kn::T1 = 0.05
    Ma::T1 = 0.0
    Pr::T1 = 2 / 3
    K::T1 = 1.0
    γ::T1 = 5 / 3
    ω::T1 = 0.81
    αᵣ::T1 = 1.0
    ωᵣ::T1 = 0.5
    μᵣ::T1 = ref_vhs_vis(Kn, αᵣ, ωᵣ)
    CFL::T1 = 0.5
    Δt::T1 = 0.01
    sim_time::T1 = 0.0
end
function Global_Data()
    geometry = [0.0, 1.0, 0.0, 1.0] # x_min,x_max,y_min,y_max
    trees_num = [20, 20]
    quadrature = [-5.0, 5.0, -5.0, 5.0]
    vs_trees_num = [20, 20] # should be even
    gas = Gas()
    ic = [1.0, 0.0, 0.0, 1.0] # prim_ini
    bc = hcat(
        [1.0, 0.0, 0.0, 1.0],
        [1.0, 0.0, 0.0, 1.0],
        [1.0, 0.0, 0.0, 1.0],
        [1.0, 0.8, 0.0, 1.0],
    ) # prim_boundary
    return Global_Data(
        geometry,
        trees_num,
        gas,
        ic,
        bc,
        quadrature,
        vs_trees_num,
        Ptr{p4est_ghost_t}(C_NULL),
        Ptr{p4est_mesh_t}(C_NULL),
        0,
        0.0,
        0,
        0,
        0,
        0,
        0,
    )
end
mutable struct Reconstruct
    pdp::PointerWrapper{PS_Data}
    dir::Integer # 1:x, 2:y, 3:z
end
mutable struct Ghost_Exchange{T1,T2,T3,T4,T5,T6}
    ghost_datas::T1
    ghost_slopes::T2
    ghost_structures::T3
    mirror_data_pointers::T4
    mirror_slope_pointers::T5
    mirror_structure_pointers::T6
end
mutable struct Trees
    data::Vector{Vector{PS_Data}}
    offset::Int
end
struct Data_Set
    PS_time_interval::Float64
    VS_time_interval::Float64
    end_time::Float64
    VS_end_time::Float64
end
function Data_Set()
    return Data_Set(0.1, 0.05, 0.02, 0.1)
end
mutable struct DVM_Data
    global_data::Global_Data
    ghost_exchange::Ghost_Exchange
    ghost_wrap::Array{Ghost_PS_Data}
    trees::Trees
    faces::Array{Face}
    p4est::Ptr{p4est_t}
    data_set::Data_Set
end
mutable struct Transfer_Data{T1,T2,T3,T4}
    w::T1
    vs_levels::T2
    vs_midpoints::T3
    vs_df::T4
    # vs_sdf::T5
end
function Transfer_Data(ps_num::Int, total_vs_num::Int)
    w = Vector{Float64}(undef, (DIM + 2) * ps_num)
    vs_levels = Vector{Int8}(undef, total_vs_num)
    vs_midpoints = Vector{Float64}(undef, DIM * total_vs_num)
    vs_df = Vector{Float64}(undef, NDF * total_vs_num)
    # vs_sdf = Vector{Float64}(undef,NDF*DIM*total_vs_num)
    return Transfer_Data(w, vs_levels, vs_midpoints, vs_df)
end
# function Transfer_Data_T(ps_num::Int,total_vs_num::Int)
#     return Transfer_Data{NTuple{(DIM+2)*ps_num,Cdouble},NTuple{total_vs_num,Int8},
#         NTuple{DIM*total_vs_num,Cdouble},NTuple{NDF*total_vs_num,Cdouble},NTuple{NDF*DIM*total_vs_num,Cdouble}}
# end
# function Transfer_Data_T(ps_num::Int,total_vs_num::Int)
#     return Transfer_Data{SVector{(DIM+2)*ps_num,Cdouble},SVector{total_vs_num,Int8},
#     SVector{DIM*total_vs_num,Cdouble},SVector{NDF*total_vs_num,Cdouble},SVector{NDF*DIM*total_vs_num,Cdouble}}
# end
mutable struct Transfer_Init
    up_num::Int
    down_num::Int
    up_data::Vector
    down_data::Vector
    old_flt::Cint # old first_local_tree
    old_llt::Cint # old last_local_tree
    up_index::Int
    up_insert_index::Int
    down_index::Int
    DVM_data::DVM_Data
end

struct Solution
    midpoint::SVector{DIM,Cdouble}
    w::SVector{DIM + 2,Cdouble}
    prim::SVector{DIM + 2,Cdouble}
    qf::SVector{DIM,Cdouble}
end
mutable struct Solution_Collect{T1,T2}
    solutions::T1
    index::T2
end
struct Mesh_Points{T}
    x::T
    y::T
end
struct Field
    midpoint::Vector{Float64}
    w::Vector{Float64}
    prim::Vector{Float64}
    qf::Vector{Float64}
end
struct PS_Result
    sim_time::Float64
    field::Vector{Field}
    mesh::Mesh_Points
end
struct VS_Result
    sim_time::Float64
    ps_midpoint::Vector{Float64}
    ps_ds::Vector{Float64}
    level::Vector{Int}
    midpoint::Matrix{Float64}
    df::Matrix{Float64}
end
struct Result_Set
    global_data::Global_Data
    PS_time_interval::Float64
    VS_time_interval::Float64
    end_time::Float64
    VS_end_time::Float64
    mpi_size::Int
end
