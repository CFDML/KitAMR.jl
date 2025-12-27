struct Uniform<:AbstractInitCondType 
    ic::AbstractVector
end
struct PCoordFn<:AbstractInitCondType # A function that accepts physical coordinates and returns initial primary variable at the position.
    PCIC_fn::Function
end


struct UGKS<:AbstractDVMFluxType end
struct DVM<:AbstractDVMFluxType end
struct CAIDVM<:AbstractDVMFluxType end # Conserved DVM, which encures energy conserved flux.
const MicroFlux = Union{DVM}
const HybridFlux = Union{UGKS,CAIDVM} # to add more...


struct Rungekuta{O} <: AbstractTimeMarchingType end
struct Euler <:AbstractTimeMarchingType end
struct UGKS_Marching<:AbstractTimeMarchingType end
struct CAIDVM_Marching<:AbstractTimeMarchingType end


struct Solver
    CFL::Float64
    AMR_PS_MAXLEVEL::Int
    AMR_DYNAMIC_PS_MAXLEVEL::Int
    AMR_VS_MAXLEVEL::Int
    flux::AbstractFluxType
    time_marching::AbstractTimeMarchingType
    PS_DYNAMIC_AMR::Bool
    VS_DYNAMIC_AMR::Bool
    AMR_PS_SMOOTH::Bool
    AMR_VS_SMOOTH::Bool
end
function Solver(config::Dict)
    return Solver(config[:CFL],config[:AMR_PS_MAXLEVEL],
        haskey(config,:AMR_DYNAMIC_PS_MAXLEVEL) ? config[:AMR_DYNAMIC_PS_MAXLEVEL] : config[:AMR_PS_MAXLEVEL],
        config[:AMR_VS_MAXLEVEL],config[:flux],config[:time_marching],
        (haskey(config,:PS_DYNAMIC_AMR) ? config[:PS_DYNAMIC_AMR] : true),
        (haskey(config,:VS_DYNAMIC_AMR) ? config[:VS_DYNAMIC_AMR] : true),
        (haskey(config,:AMR_PS_SMOOTH) ? config[:AMR_PS_SMOOTH] : true),
        (haskey(config,:AMR_VS_SMOOTH) ? config[:AMR_VS_SMOOTH] : true)
        )
end
mutable struct UDF
    static_ps_refine_flag::Function
    dynamic_ps_refine_flag::Function
    static_ps_coarsen_flag::Function
    vs_refine_flag::Function
    vs_coarsen_flag::Function 
    UDF()=new()
end
null_udf(args...;kwargs...) = false
mutable struct Output
    vtk_celltype::DataType
    vs_vtk_celltype::DataType
    anim_dt::Float64
    vs_output_criterion::Function
    anim_index::Int
    Output(config::Dict) = (n = new(config[:DIM]==2 ? Triangle : Tetra,
    config[:DIM]==2 ? Pixel : Voxel,
    0.,null_udf,0
    );
        for i in fieldnames(Output)
            if haskey(config,i)
                setfield!(n,i,config[i])
            end
        end;
    n
    )
end
struct Configure{DIM,NDF}
    geometry::Vector{Float64}
    trees_num::Vector{Int64}
    quadrature::Union{AbstractQuadrature,Vector{Float64}}
    vs_trees_num::Vector{Int64}
    IC::AbstractInitCondType
    domain::Vector{Domain}
    IB::Vector{AbstractBoundary}
    gas::Gas
    solver::Solver
    output::Output
    user_defined::UDF
end
function config_IB(ib::Circle,config::Dict)
    ds = minimum([(config[:geometry][2i]-config[:geometry][2i-1])/config[:trees_num][i]/2^config[:AMR_PS_MAXLEVEL] for i in 1:config[:DIM]])
    Circle(ib,ds)
end
function config_IB(ib::Sphere,config::Dict)
    ds = minimum([(config[:geometry][2i]-config[:geometry][2i-1])/config[:trees_num][i]/2^config[:AMR_PS_MAXLEVEL] for i in 1:config[:DIM]])
    Sphere(ib,ds)
end
function config_IB(ib::Vertices,config::Dict)
    Vertices(ib,config)
end
function Configure(config::Dict)
    gas = Gas()
    for i in fieldnames(Gas)
        if haskey(config,i)
            setfield!(gas,i,config[i])
        end
    end
    if !haskey(config,:μᵣ)
        gas.μᵣ = ref_vhs_vis(gas.Kn,gas.αᵣ,gas.ωᵣ)
    end
    IB = haskey(config,:IB) ? config[:IB] : []
    for i in eachindex(IB)
        IB[i] = config_IB(IB[i],config)
    end
    user_defined = UDF()
    for i in fieldnames(UDF)
        if haskey(config,i)
            setfield!(user_defined,i,config[i])
        else
            setfield!(user_defined,i,null_udf)
        end
    end
    return Configure{config[:DIM],config[:NDF]}(config[:geometry],config[:trees_num],
        config[:quadrature],config[:vs_trees_num],config[:IC],config[:domain],IB,
        gas,Solver(config),Output(config),user_defined)
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
mutable struct Residual
    step::Int
    residual::Vector{Float64}
    sumRes::Vector{Float64}
    sumAvg::Vector{Float64}
    redundant_step::Int
end
function Residual(DIM::Int)
    return Residual(1,ones(DIM+2),zeros(DIM+2),zeros(DIM+2),0)
end
mutable struct Status
    max_vs_num::Int # maximum vs_num among ghost quadrants
    gradmax::Vector{Float64}
    Δt::Float64
    Δt_ξ::Float64
    sim_time::Float64
    ps_adapt_step::Int
    vs_adapt_step::Int
    partition_step::Int
    residual::Residual
    save_flag::Base.RefValue{Bool}
    ib_ghost::Bool
end
function Status(config)
    DIM = config[:DIM]
    trees_num = config[:trees_num]
    geometry = config[:geometry]
    vs_trees_num = config[:vs_trees_num]
    quadrature = config[:quadrature]
    ds = [(geometry[2*i]-geometry[2*i-1])/trees_num[i]/2^config[:AMR_PS_MAXLEVEL] for i in 1:DIM]
    U = isa(quadrature,Vector) ? [max(quadrature[2*i],abs(quadrature[2*i-1])) -
        (quadrature[2*i] - quadrature[2*i-1]) / vs_trees_num[i]/
        2^config[:AMR_VS_MAXLEVEL] / 2 for i in 1:DIM] : [maximum(abs.(quadrature.vcoords)) for _ in 1:DIM]
    Δt = config[:CFL]*minimum(ds ./ U)
    return Status(0,ones(DIM+2),Δt*TIME_STEP_CONTRACT_RATIO,Δt,0.,1,1,1,Residual(DIM),Ref(false),false)
end

mutable struct Global_Data{DIM,NDF}
    config::Configure{DIM,NDF}
    forest::Forest{DIM}
    status::Status
    Global_Data(config::Dict) = (n = new{config[:DIM],config[:NDF]}();
    n.config = Configure(config);
    n.forest = Forest(config[:DIM]);
    n.status = Status(config);
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

mutable struct Field{DIM,NDF}
    trees::PS_Trees{DIM,NDF}
    faces::Vector{AbstractFace}
    immersed_boundary::ImmersedBoundary
end


mutable struct AMR{DIM,NDF}
    global_data::Global_Data
    ghost::Ghost
    field::Field{DIM,NDF}
    AMR(global_data::Global_Data{DIM,NDF},field) where{DIM,NDF} = (n = new{DIM,NDF}();
        n.global_data = global_data;
        n.field = field;
        n
    )
end
# partition

struct Transfer_Data{DIM,NDF}
    encs::Vector{Int}
    w::Vector{Float64}
    vs_levels::Vector{Int8}
    vs_midpoints::Vector{Float64}
    vs_df::Vector{Float64}
end
function Transfer_Data(DIM::Integer,NDF::Integer,ps_num::Integer,total_vs_num::Integer)
    return Transfer_Data{DIM,NDF}(
        Vector{Int}(undef, (SOLID_CELL_ID_NUM+1)*ps_num),
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
