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
    domain::Vector{Domain}
    IB::Vector{AbstractBoundary}
    IB_sort::Symbol
    IB_interp::Symbol
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
    domain = Domain[]
    bc = AbstractBoundary[]
    for i = 1:2*config[:DIM]
        push!(domain,Domain{config[:domaintype][i]}(i,config[:domainbc][i]))
    end
    if !isempty(config[:boundarydefine])
        for i in eachindex(config[:boundarydefine])
            push!(bc,config[:boundarydefine][i])
        end
    end
    return Configure{config[:DIM],config[:NDF]}(config[:geometry],config[:trees_num],
        config[:quadrature],config[:vs_trees_num],config[:ic],domain,bc,config[:ib_sort],
        config[:ib_interp],gas,Solver(config))
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
    data::Vector{Vector{AbstractPsData{DIM,NDF}}}
    offset::Int
end

mutable struct Ghost
    ghost_exchange::Ghost_Exchange
    ghost_wrap::Vector{AbstractGhostPsData}
end

mutable struct Boundary{DIM,NDF}
    solid_cells::Vector{SolidCells{DIM,NDF}} # Element corresponds to one IB boundary
    Numbers::Vector{Vector{Int}}
    image_points::Vector{Vector{Vector{Float64}}} # Image points of solid_cells sorted by quadid
    aux_points::Vector{Vector{Vector{Float64}}} # Midpoints of aux_points sorted by quadid
    IB_cells::Vector{IBCells} # Element corresponds to one IB boundary
    IB_ranks_table::Vector{Vector{PS_Data{DIM,NDF}}} # Local IB nodes belonging to solidcells in different processors
    IB_buffer::IBBuffer # Buffer for IB communication. Store pointers for memory free.
end
mutable struct Field{DIM,NDF}
    trees::PS_Trees{DIM,NDF}
    faces::Vector{Face}
    boundary::Boundary{DIM,NDF}
    # solid_cells::Vector{Vector{PS_Data{DIM,NDF}}} # Outer vector corresbonds to boundaries
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
    encs::Vector{Int8}
    w::Vector{Float64}
    vs_levels::Vector{Int8}
    vs_midpoints::Vector{Float64}
    vs_df::Vector{Float64}
end
function Transfer_Data(DIM::Integer,NDF::Integer,ps_num::Integer,total_vs_num::Integer)
    return Transfer_Data{DIM,NDF}(
        Vector{Int8}(undef, 2*ps_num),
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
