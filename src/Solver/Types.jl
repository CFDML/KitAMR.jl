"""
$(TYPEDEF)

Uniform initial condition. The field will be initialized with the Maxwellian distribution determined by the unified primary macroscopic variables in ic.

## Fields

$(TYPEDFIELDS)
"""
struct Uniform<:AbstractInitCondType
    ic::AbstractVector
end
"""
$(TYPEDEF)

Physical-coordinates-determined initial condition. The field will be initialized with the Maxwellian distribution determined by the user-defined-function in PCIC_fn.
The function will return primary macroscopic variables vector according to the information of the physical grid (mostly the physical coordinates).

## Fields

$(TYPEDFIELDS)
"""
struct PCoordFn<:AbstractInitCondType # A function that accepts physical coordinates and returns initial primary variable at the position.
    PCIC_fn::Function
end



"""
$(TYPEDEF)

Structure of solver configuration.

## Fields

$(TYPEDFIELDS)
"""
struct Solver{DIM,NDF}
    """
    Courant-Friedrichs-Lewy number. Default is `0.4`.
    """
    CFL::Float64
    """
    Maximum level of the (static) AMR in physical space. **Mandatory**.
    """
    AMR_PS_MAXLEVEL::Int
    """
    Maximum level of the dynamic AMR in physical space. In most cases, it should be smaller than the static one. Default is equal to `AMR_PS_MAXLEVEL`.
    """
    AMR_DYNAMIC_PS_MAXLEVEL::Int
    """
    Maximum level of the AMR in velocity space. **Mandatory**.
    """
    AMR_VS_MAXLEVEL::Int
    """
    The numerical flux type. **Mandatory**.
    """
    flux::Type{Tf} where {Tf<:AbstractFluxType}
    """
    The time-marching scheme. **Mandatory**.
    """
    time_marching::Type{Tt} where {Tt<:AbstractTimeMarchingType}
    """
    The dynamic AMR in physical space is open or not. Default is `true`.
    """
    PS_DYNAMIC_AMR::Bool
    """
    The dynamic AMR in velocity space is open or not. Default is `true`.
    """
    VS_DYNAMIC_AMR::Bool
end
function Solver(config::Dict)
    return Solver{config[:DIM],config[:NDF]}(
        config[:CFL],
        config[:AMR_PS_MAXLEVEL],
        haskey(config, :AMR_DYNAMIC_PS_MAXLEVEL) ? config[:AMR_DYNAMIC_PS_MAXLEVEL] :
        config[:AMR_PS_MAXLEVEL],
        config[:AMR_VS_MAXLEVEL],
        config[:flux],
        config[:time_marching],
        (haskey(config, :PS_DYNAMIC_AMR) ? config[:PS_DYNAMIC_AMR] : true),
        (haskey(config, :VS_DYNAMIC_AMR) ? config[:VS_DYNAMIC_AMR] : true),
    )
end
function Solver(; kwargs...)
    CFL = haskey(kwargs, :CFL) ? kwargs[:CFL] : 0.4
    AMR_DYNAMIC_PS_MAXLEVEL =
        haskey(kwargs, :AMR_DYNAMIC_PS_MAXLEVEL) ? kwargs[:AMR_DYNAMIC_PS_MAXLEVEL] :
        kwargs[:AMR_PS_MAXLEVEL]
    PS_DYNAMIC_AMR = haskey(kwargs, :PS_DYNAMIC_AMR) ? kwargs[:PS_DYNAMIC_AMR] : true
    VS_DYNAMIC_AMR = haskey(kwargs, :VS_DYNAMIC_AMR) ? kwargs[:VS_DYNAMIC_AMR] : true
    return Solver{kwargs[:DIM],kwargs[:NDF]}(
        CFL,
        kwargs[:AMR_PS_MAXLEVEL],
        AMR_DYNAMIC_PS_MAXLEVEL,
        kwargs[:AMR_VS_MAXLEVEL],
        kwargs[:flux],
        kwargs[:time_marching],
        PS_DYNAMIC_AMR,
        VS_DYNAMIC_AMR,
    )
end

"""
$(TYPEDEF)

Structure of user-defined-functions.

## Fields

$(TYPEDFIELDS)
"""
mutable struct UDF
    """
    The static AMR flag in physical space.
    """
    static_ps_refine_flag::Function
    """
    The dynamic AMR flag in physical space.
    """
    dynamic_ps_refine_flag::Function
    """
    The static AMR flag in velocity space.
    """
    static_vs_refine_flag::Function
end
null_udf(args...; kwargs...) = false
function UDF(; kwargs...)
    static_ps_refine_flag =
        haskey(kwargs, :static_ps_refine_flag) ? kwargs[:static_ps_refine_flag] : null_udf
    dynamic_ps_refine_flag =
        haskey(kwargs, :dynamic_ps_refine_flag) ? kwargs[:dynamic_ps_refine_flag] : null_udf
    static_vs_refine_flag =
        haskey(kwargs, :static_vs_refine_flag) ? kwargs[:static_vs_refine_flag] : null_udf
    return UDF(static_ps_refine_flag, dynamic_ps_refine_flag, static_vs_refine_flag)
end
"""
$(TYPEDEF)

Structure of output information.

## Fields

$(TYPEDFIELDS)
"""
mutable struct Output
    """
    VTK cell type in physical space. Available options: `Pixel` and `Triangle` for 2D; `Voxel` and `Tetra` for 3D. Default is [`Triangle`](@ref) for 2D and [`Tetra`](@ref) for 3D.
    """
    vtk_celltype::Type{Tp} where {Tp<:AbstractVTKCellType}
    """
    VTK cell type in velocity space. Available options: `Pixel` for 2D.
    """
    vs_vtk_celltype::Type{Tv} where {Tv<:AbstractVTKCellType}
    """
    Time interval between animation frames.
    """
    anim_dt::Float64
    """
    User-defined-function determining the physical cells need to output the velocity space.
    """
    vs_output_criterion::Function
    """
    Index of the last saved animation frame.
    """
    anim_index::Int
end
function Output(config::Dict)
    output = Output(
        config[:DIM]==2 ? Triangle : Tetra,
        config[:DIM]==2 ? Pixel : Voxel,
        0.0,
        null_udf,
        0,
    )
    for i in fieldnames(Output)
        if haskey(config, i)
            setfield!(n, i, config[i])
        end
    end
    return output
end
function Output(solver::Solver{DIM,NDF}; kwargs...) where {DIM,NDF}
    output = Output(DIM==2 ? Triangle : Tetra, DIM==2 ? Pixel : Voxel, 0.0, null_udf, 0)
    for i in fieldnames(Output)
        if haskey(kwargs, i)
            setfield!(n, i, kwargs[i])
        end
    end
    return output
end

"""
$(TYPEDEF)

Structure of problem configuration. `DIM` and `NDF` represent the dimension number of the problem and the number of the distribution function.
This struct plays an key role in the solution process.

## Fields

$(TYPEDFIELDS)
"""
struct Configure{DIM,NDF}
    """
    Range of the simulated domain. As an example, for 2D case, it should be aligned as [xmin,xmax,ymin,ymax].
    """
    geometry::Vector{Float64}
    """
    Number of tree roots for each dimension. For 2D case, it should be aligned as [x_num,y_num]
    """
    trees_num::Vector{Int64}
    """
    For Vector type, it represents the range of the velocity space. For 2D case, it should be aligned as [umin,umax,vmin,vmax]. The `AbstractQuadrature` type is reserved for other quadrature rule like Gauss-Hermite.
    """
    quadrature::Union{AbstractQuadrature,Vector{Float64}}
    """
    Number of tree roots for each dimension in velocity space.
    """
    vs_trees_num::Vector{Int64}
    """
    Initial condition defined by [`AbstractInitCondType`](@ref).
    """
    IC::AbstractInitCondType
    """
    Types of the domain boundary. For 2D case, the vector should catain 4 elements corresponding to the 4 domain boundaries. The element type is defined by [`Domain`](@ref).
    """
    domain::Vector{Domain}
    """
    Immersed boundary. Multiple boundaries are supported. The element type is defined by [`AbstractBoundaryType`](@ref)
    """
    IB::Vector{AbstractBoundaryType}
    """
    Properties of the simulated gas defined by [`Gas`](@ref).
    """
    gas::Gas
    """
    Setup of the solver defined by [`Solver`](@ref).
    """
    solver::Solver{DIM,NDF}
    """
    Setup of the output form defined by [`Output`](@ref).
    """
    output::Output
    """
    Functions defined by users, including some criteria.
    """
    user_defined::UDF
end
function config_IB(ib::Circle, config::Dict)
    ds = minimum([
        (config[:geometry][2i]-config[:geometry][2i-1])/config[:trees_num][i]/2^config[:AMR_PS_MAXLEVEL]
        for i = 1:config[:DIM]
    ])
    Circle(ib, ds)
end
function config_IB(ib::Sphere, config::Dict)
    ds = minimum([
        (config[:geometry][2i]-config[:geometry][2i-1])/config[:trees_num][i]/2^config[:AMR_PS_MAXLEVEL]
        for i = 1:config[:DIM]
    ])
    Sphere(ib, ds)
end
function config_IB(ib::Vertices, config::Dict)
    Vertices(ib, config)
end
function config_IB(ib::Triangles, config::Dict)
    Triangles(ib, config)
end
function Configure(config::Dict)
    gas = Gas()
    for i in fieldnames(Gas)
        if haskey(config, i)
            setfield!(gas, i, config[i])
        end
    end
    if !haskey(config, :μᵣ)
        gas.μᵣ = ref_vhs_vis(gas.Kn, gas.αᵣ, gas.ωᵣ)
    end
    IB = haskey(config, :IB) ? config[:IB] : []
    for i in eachindex(IB)
        IB[i] = config_IB(IB[i], config)
    end
    user_defined = UDF()
    for i in fieldnames(UDF)
        if haskey(config, i)
            setfield!(user_defined, i, config[i])
        else
            setfield!(user_defined, i, null_udf)
        end
    end
    return Configure{config[:DIM],config[:NDF]}(
        config[:geometry],
        config[:trees_num],
        config[:quadrature],
        config[:vs_trees_num],
        config[:IC],
        config[:domain],
        IB,
        gas,
        Solver(config),
        Output(config),
        user_defined,
    )
end
function Configure(solver::Solver{DIM,NDF}; kwargs...) where {DIM,NDF}
    geometry = kwargs[:geometry]
    trees_num = kwargs[:trees_num]
    AMR_PS_MAXLEVEL = solver.AMR_PS_MAXLEVEL
    IB = haskey(kwargs, :IB) ? kwargs[:IB] : []
    config_dict = Dict(
        :DIM=>DIM,
        :geometry=>geometry,
        :trees_num=>trees_num,
        :AMR_PS_MAXLEVEL=>AMR_PS_MAXLEVEL,
    )
    for i in eachindex(IB)
        IB[i] = config_IB(IB[i], config_dict)
    end
    return Configure{DIM,NDF}(
        geometry,
        trees_num,
        kwargs[:quadrature],
        kwargs[:vs_trees_num],
        kwargs[:IC],
        kwargs[:domain],
        IB,
        kwargs[:gas],
        solver,
        kwargs[:output],
        kwargs[:user_defined],
    )
end
"""
$(TYPEDEF)
Structure related to `p4est`.
$(TYPEDFIELDS)
"""
mutable struct Forest{DIM}
    p4est::P_pxest_t
    ghost::P_pxest_ghost_t
    mesh::P_pxest_mesh_t
    Forest(DIM) = (
        n = new{DIM}();
        n.p4est = Ptr{pxest_ts[DIM-1]}(C_NULL);
        n.ghost = Ptr{pxest_ghost_ts[DIM-1]}(C_NULL);
        n.mesh = Ptr{pxest_mesh_ts[DIM-1]}(C_NULL);
        n
    )
end
"""
$(TYPEDEF)
$(TYPEDFIELDS)
"""
mutable struct Residual
    """
    Current step number.
    """
    step::Int
    """
    Residual of primary variables, is updated by [`residual_check!`](@ref).
    """
    residual::Vector{Float64}
    """
    Summation of residuals.
    """
    sumRes::Vector{Float64}
    """
    Summation of primary variables to average `sumRes`.
    """
    sumAvg::Vector{Float64}
    """
    Count of redundant step to insure the convergence.
    """
    redundant_step::Int
end
function Residual(DIM::Int)
    return Residual(1, ones(DIM+2), zeros(DIM+2), zeros(DIM+2), 0)
end

"""
$(TYPEDEF)

Structure of simulation status. It contains the real time information, and is updated every step.

## Fields

$(TYPEDFIELDS)
"""
mutable struct Status
    """
    Maximum number of the velocity grids among ghost quadrants. # maximum vs_num among ghost quadrants
    """
    max_vs_num::Int # maximum vs_num among ghost quadrants
    """
    Maximum absolute value of the gradients of conserved variables.
    """
    gradmax::Vector{Float64}
    """
    Maximum value of conserved variables.
    """
    wmax::Vector{Float64}
    """
    Minimum value of conserved variables.
    """
    wmin::Vector{Float64}
    """
    Time step size used in [`iterate!`](@ref).
    """
    Δt::Float64
    """
    Time step size constraint by grid size.
    """
    Δt_ξ::Float64
    """
    Dimensionless simulation time.
    """
    sim_time::Float64
    """
    Number of steps after last AMR in physical space.
    """
    ps_adapt_step::Int
    """
    Number of steps after last AMR in velocity space.
    """
    vs_adapt_step::Int
    """
    Number of steps after last partition.
    """
    partition_step::Int
    """
    Residual of conserved variables defined by [`Residual`](@ref).
    """
    residual::Residual
    """
    Flag indicating whether to save.
    """
    save_flag::Base.RefValue{Bool}
end
function Status(config::Dict)
    DIM = config[:DIM]
    trees_num = config[:trees_num]
    geometry = config[:geometry]
    vs_trees_num = config[:vs_trees_num]
    quadrature = config[:quadrature]
    ds = [
        (geometry[2*i]-geometry[2*i-1])/trees_num[i]/2^config[:AMR_PS_MAXLEVEL] for
        i = 1:DIM
    ]
    U =
        isa(quadrature, Vector) ?
        [
            max(quadrature[2*i], abs(quadrature[2*i-1])) -
            (quadrature[2*i] - quadrature[2*i-1]) / vs_trees_num[i] /
            2^config[:AMR_VS_MAXLEVEL] / 2 for i = 1:DIM
        ] : [maximum(abs.(quadrature.vcoords)) for _ = 1:DIM]
    Δt_ξ = config[:CFL]*minimum(ds ./ U)
    return Status(
        0,
        zeros(DIM+2),
        zeros(DIM+2),
        zeros(DIM+2),
        Δt_ξ,
        Δt_ξ,
        0.0,
        1,
        1,
        1,
        Residual(DIM),
        Ref(false),
    )
end
function Status(config::Configure{DIM,NDF}) where {DIM,NDF}
    trees_num = config.trees_num
    geometry = config.geometry
    vs_trees_num = config.vs_trees_num
    quadrature = config.quadrature
    ds = [
        (geometry[2*i]-geometry[2*i-1])/trees_num[i]/2^config.solver.AMR_PS_MAXLEVEL for
        i = 1:DIM
    ]
    U =
        isa(quadrature, Vector) ?
        [
            max(quadrature[2*i], abs(quadrature[2*i-1])) -
            (quadrature[2*i] - quadrature[2*i-1]) / vs_trees_num[i] /
            2^config.solver.AMR_VS_MAXLEVEL / 2 for i = 1:DIM
        ] : [maximum(abs.(quadrature.vcoords)) for _ = 1:DIM]
    Δt_ξ = config.solver.CFL*minimum(ds ./ U)
    return Status(
        0,
        zeros(DIM+2),
        zeros(DIM+2),
        zeros(DIM+2),
        Δt_ξ,
        Δt_ξ,
        0.0,
        1,
        1,
        1,
        Residual(DIM),
        Ref(false),
    )
end

"""
$(TYPEDEF)
Structure of required information during a simulation.
$(TYPEDFIELDS)
"""
mutable struct Global_Data{DIM,NDF}
    """
    Defined by [`Configure`](@ref).
    """
    config::Configure{DIM,NDF}
    """
    Defined by [`Forest`](@ref).
    """
    forest::Forest{DIM}
    """
    Defined by [`Status`](@ref).
    """
    status::Status
    Global_Data(config::Dict) = (
        n = new{config[:DIM],config[:NDF]}();
        n.config = Configure(config);
        n.forest = Forest(config[:DIM]);
        n.status = Status(config);
        n
    )
    Global_Data(config::Configure{DIM,NDF}) where {DIM,NDF} = (
        n = new{DIM,NDF}();
        n.config = config;
        n.forest = Forest(DIM);
        n.status = Status(config);
        n
    )
end

mutable struct P4est_PS_Data
    ps_data::Ptr{Nothing}
end

"""
$(TYPEDEF)
Structure of ghost data and pointers used by `p4est`.
$(TYPEDFIELDS)
"""
mutable struct Ghost_Exchange
    ghost_datas::Any
    ghost_slopes::Any
    ghost_structures::Any
    mirror_data_pointers::Any
    mirror_slope_pointers::Any
    mirror_structure_pointers::Any
end

"""
$(TYPEDEF)
Structure of cells data managed by `Julia`.
$(TYPEDFIELDS)
"""
mutable struct PS_Trees{DIM,NDF}
    data::Vector{Vector{AbstractPsData{DIM,NDF}}}
    """
    Offset of treeid used in partition.
    """
    offset::Int
end

"""
$(TYPEDEF)
Structure of ghost layer for communication between processors.
$(TYPEDFIELDS)
"""
mutable struct Ghost
    ghost_exchange::Ghost_Exchange
    ghost_wrap::Vector{AbstractGhostPsData}
end

"""
$(TYPEDEF)
Structure of fields data.
$(TYPEDFIELDS)
"""
mutable struct Field{DIM,NDF}
    """
    Cells data defined by [`PS_Trees`](@ref).
    """
    trees::PS_Trees{DIM,NDF}
    """
    Mapping between faces and cells. The element type is defined by [`AbstractFace`](@ref).
    """
    faces::Vector{AbstractFace}
    immersed_boundary::ImmersedBoundary
end

"""
$(TYPEDEF)
Structure that collects all of the data.
$(TYPEDFIELDS)
"""
mutable struct KitAMR_Data{DIM,NDF}
    global_data::Global_Data
    ghost::Ghost
    field::Field{DIM,NDF}
    KitAMR_Data(global_data::Global_Data{DIM,NDF}, field) where {DIM,NDF} =
        (n = new{DIM,NDF}(); n.global_data = global_data; n.field = field; n)
end
# partition

struct Transfer_Data{DIM,NDF}
    encs::Vector{Int}
    w::Vector{Float64}
    vs_levels::Vector{Int8}
    vs_midpoints::Vector{Float64}
    vs_df::Vector{Float64}
end
function Transfer_Data(DIM::Integer, NDF::Integer, ps_num::Integer, total_vs_num::Integer)
    return Transfer_Data{DIM,NDF}(
        Vector{Int}(undef, (SOLID_CELL_ID_NUM+1)*ps_num),
        Vector{Float64}(undef, (DIM+2) * ps_num),
        Vector{Int8}(undef, total_vs_num),
        Vector{Float64}(undef, DIM * total_vs_num),
        Vector{Float64}(undef, NDF * total_vs_num),
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
end
