pxest_t = Union{p4est_t,p8est_t}
pxest_ghost_t = Union{p4est_ghost_t,p8est_ghost_t}
pxest_mesh_t = Union{p4est_mesh_t,p8est_mesh_t}

P_pxest_t = Union{Ptr{p4est_t}, Ptr{p8est_t}}
P_pxest_ghost_t = Union{Ptr{p4est_ghost_t}, Ptr{p8est_ghost_t}}
P_pxest_mesh_t = Union{Ptr{p4est_mesh_t}, Ptr{p8est_mesh_t}}
P_pxest_quadrant_t = Union{Ptr{p4est_quadrant_t}, Ptr{p8est_quadrant_t}}
P_pxest_iter_volume_info_t = Union{Ptr{p4est_iter_volume_info_t}, Ptr{p8est_iter_volume_info_t}}
P_pxest_iter_face_info_t = Union{Ptr{p4est_iter_face_info_t}, Ptr{p8est_iter_face_info_t}}

PW_pxest_t = Union{PointerWrapper{p4est_t}, PointerWrapper{p8est_t}}
PW_pxest_ghost_t = Union{PointerWrapper{p4est_ghost_t}, PointerWrapper{p8est_ghost_t}}
PW_pxest_mesh_t = Union{PointerWrapper{p4est_mesh_t}, PointerWrapper{p8est_mesh_t}}
PW_pxest_quadrant_t = Union{PointerWrapper{p4est_quadrant_t}, PointerWrapper{p8est_quadrant_t}}
PW_pxest_iter_volume_info_t = Union{PointerWrapper{p4est_iter_volume_info_t}, PointerWrapper{p8est_iter_volume_info_t}}
PW_pxest_iter_face_info_t = Union{PointerWrapper{p4est_iter_face_info_t},PointerWrapper{p8est_iter_face_info_t}}

struct Solver
    CFL::Float64
    AMR_PS_MAXLEVEL::Int
    AMR_VS_MAXLEVEL::Int
end
function Solver(config::Dict)
    return Solver(config[:CFL],config[:AMR_PS_MAXLEVEL],config[:AMR_VS_MAXLEVEL])
end

struct Configure{DIM,NDF}
    geometry::Vector{Float64}
    trees_num::Vector{Int64}
    quadrature::Vector{Float64}
    vs_trees_num::Vector{Int64}
    ic::Vector{Float64}
    bc::Matrix{Float64}
    gas::Gas
    solver::Solver
end
function Configure(config::Dict)
    gas = Gas()
    for i in fieldnames(Gas)
        if haskey(config,i)
            setfield!(gas,i,config[i])
        end
    end;
    return Configure{config[:DIM],config[:NDF]}(config[:geometry],config[:trees_num],config[:quadrature],config[:vs_trees_num],config[:ic],config[:bc],gas,Solver(config))
end

mutable struct Forest{DIM}
    p4est::P_pxest_t
    ghost::P_pxest_ghost_t
    mesh::P_pxest_mesh_t
    Forest(DIM) = (n = new{DIM}();
    n.p4est = Ptr{pxest_ts[DIM-1]}(C_NULL);
    n.ghost = Ptr{pxest_ghost_ts[DIM-1]}(C_NULL);
    n.mesh = Ptr{pxest_mesh_ts[DIM-1]}(C_NULL);
    n
    )
end

mutable struct Status
    max_vs_num::Int
    gradmax::Float64
    Δt::Float64
    sim_time::Float64
    ps_adapt_step::Int
    vs_adapt_step::Int
    partition_step::Int
end
function Status()
    return Status(0,0.,1.,0.,0,0,0)
end

mutable struct Global_Data{DIM,NDF}
    config::Configure{DIM,NDF}
    forest::Forest{DIM}
    status::Status
    Global_Data(config::Dict) = (n = new{config[:DIM],config[:NDF]}();
    n.config = Configure(config);
    n.forest = Forest(config[:DIM]);
    n.status = Status();
    n
    )
end

mutable struct P4est_PS_Data
    ps_data::Ptr{Nothing}
end

mutable struct Ghost_Exchange
    ghost_datas
    ghost_slopes
    ghost_structures
    mirror_data_pointers
    mirror_slope_pointers
    mirror_structure_pointers
end

mutable struct PS_Trees{DIM,NDF} 
    data::Vector{Vector{PS_Data{DIM,NDF}}}
    offset::Int
end

mutable struct Ghost
    ghost_exchange::Ghost_Exchange
    ghost_wrap::Vector{Ghost_PS_Data}
end

mutable struct Field{DIM,NDF}
    trees::PS_Trees{DIM,NDF}
    faces::Vector{Face}
end

mutable struct AMR{DIM,NDF}
    global_data::Global_Data
    ghost::Ghost
    field::Field{DIM,NDF}
end
function AMR(global_data::Global_Data{DIM,NDF},ghost::Ghost,field::Field{DIM,NDF}) where{DIM,NDF}
    return AMR{DIM,NDF}(global_data,ghost,field)
end

# partition

mutable struct Transfer_Data{DIM,NDF}
    w::Vector{Float64}
    vs_levels::Vector{Int8}
    vs_midpoints::Vector{Float64}
    vs_df::Vector{Float64}
end
function Transfer_Data(DIM::Integer,NDF::Integer,ps_num::Integer,total_vs_num::Integer)
    return Transfer_Data{DIM,NDF}(
        Vector{Float64}(undef, (DIM+2) * ps_num),
        Vector{Int8}(undef, total_vs_num),
        Vector{Float64}(undef, DIM * total_vs_num),
        Vector{Float64}(undef, NDF * total_vs_num)
    )
end

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
    amr::AMR
end

struct Solution_2D
    midpoint::SVector{2,Float64}
    w::SVector{4,Float64}
    prim::SVector{4,Float64}
    qf::SVector{2,Float64}
end
struct Solution_3D
    midpoint::SVector{3,Float64}
    w::SVector{5,Float64}
    prim::SVector{5,Float64}
    qf::SVector{3,Float64}
end

struct SolverSet
    MPI_size::Int
    global_data::Global_Data
end