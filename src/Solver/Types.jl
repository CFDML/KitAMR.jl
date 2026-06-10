"""
$(TYPEDEF)

Uniform initial condition: every cell is initialized to the Maxwellian distribution determined
by one primitive-variable vector `ic`, identical across the whole domain.

The primitive vector is ordered `[ρ, U₁, …, U_DIM, λ]` (length `DIM+2`), where `ρ` is the
density, `Uᵢ` the bulk-velocity components, and `λ = ρ/(2p)` the inverse temperature
(`p` is the pressure). This ordering is shared by [`PCoordFn`](@ref) and the [`Domain`](@ref)
boundary states.

## Fields

$(TYPEDFIELDS)
"""
struct Uniform<:AbstractInitCond
    """
    Primitive macroscopic variables `[ρ, U₁, …, U_DIM, λ]` applied uniformly to every cell.
    """
    ic::AbstractVector
end
"""
$(TYPEDEF)

Coordinate-dependent initial condition: each cell is initialized to the Maxwellian determined
by the primitive vector returned by the user function `PCIC_fn` evaluated at the cell centre.
Use this for any non-uniform initial state (shock tubes, Riemann/blast initial data, buffer
zones around immersed bodies, …).

`PCIC_fn` must have the signature

    PCIC_fn(midpoint::Vector{Float64}, kinfo::KInfo) -> Vector{Float64}

where `midpoint` is the cell-centre coordinate (length `DIM`) and the return value is the
primitive vector `[ρ, U₁, …, U_DIM, λ]` (length `DIM+2`, `λ = ρ/(2p)`). See the
[User-defined functions](@ref) page for the full convention and worked examples.

## Fields

$(TYPEDFIELDS)
"""
struct PCoordFn<:AbstractInitCond # A function that accepts physical coordinates and returns initial primary variable at the position.
    """
    Function `(midpoint, kinfo) -> prim` returning the primitive vector at a physical coordinate.
    """
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
    Numerical flux type. **Mandatory**.
    """
    flux::Type{Tf} where {Tf<:AbstractFluxType}
    """
    Time-marching scheme. **Mandatory**.
    """
    time_marching::Type{Tt} where {Tt<:AbstractTimeMarchingType}
    """
    Dynamic AMR in physical space is open or not. Default is `true`.
    """
    PS_DYNAMIC_AMR::Bool
    """
    Dynamic AMR in velocity space is open or not. Default is `true`.
    """
    VS_DYNAMIC_AMR::Bool
    """
    Criterion of AMR in L\"ohner criterion of physical space. Default is `0.2`.
    """
    ADAPT_COEFFI_PS::Float64
    """
    Redundancy coefficient of AMR in velocity space. Default is `0.125`.
    """
    ADAPT_COEFFI_VS_GLOBAL::Float64
    """
    Proportion that a refinement-required velocity cell contributes to the macroscopic quantities in a local physical cell. Default is `0.01`.
    """
    ADAPT_COEFFI_VS_LOCAL::Float64
    """
    Velocity-space refinement criterion. `:lohner` (default) uses the moment-weighted Löhner indicator with `local_contribution_*` as a relative mass/energy floor; `:contribution` uses the legacy magnitude/contribution flags.
    """
    ADAPT_VS_MODE::Symbol
    """
    Threshold of the moment-weighted Löhner indicator (in `[0,1]`) for velocity-space refinement when `ADAPT_VS_MODE == :lohner`. Default is `0.6`.
    """
    ADAPT_COEFFI_VS_LOHNER::Float64
    """
    Tolerance for convergence judgement. Default is `1e-6`.
    """
    TOLERANCE::Float64
    """
    Number of steps between two checks of status. Default is `100`.
    """
    ST_CHECK_INTERVAL::Int
    """
    Number of redundant steps after convergence criterion has been satisfied. Default is `100`.
    """
    REDUNDANT_STEPS_NUM::Int
    """
    Simulation termination time (dimensionless). The main loop runs until `sim_time` reaches this value, and the last time step is shrunk to land exactly on it (see [`limit_Δt!`](@ref) / [`reached_max_time`](@ref)). Default is `Inf` (no time limit; terminate by convergence).
    """
    max_sim_time::Float64
end
function Solver(config::Dict)
    ADAPT_COEFFI_PS = haskey(config, :ADAPT_COEFFI_PS) ? config[:ADAPT_COEFFI_PS] : 0.25
    ADAPT_COEFFI_VS_GLOBAL =
        haskey(config, :ADAPT_COEFFI_VS_GLOBAL) ? config[:ADAPT_COEFFI_VS_GLOBAL] : 0.125
    ADAPT_COEFFI_VS_LOCAL =
        haskey(config, :ADAPT_COEFFI_VS_LOCAL) ? config[:ADAPT_COEFFI_VS_LOCAL] : 1e-2
    ADAPT_VS_MODE =
        haskey(config, :ADAPT_VS_MODE) ? Symbol(config[:ADAPT_VS_MODE]) : :lohner
    ADAPT_COEFFI_VS_LOHNER =
        haskey(config, :ADAPT_COEFFI_VS_LOHNER) ? config[:ADAPT_COEFFI_VS_LOHNER] : 0.6
    TOLERANCE = haskey(config, :TOLERANCE) ? config[:TOLERANCE] : 1e-6
    ST_CHECK_INTERVAL =
        haskey(config, :ST_CHECK_INTERVAL) ? config[:ST_CHECK_INTERVAL] : 100
    REDUNDANT_STEPS_NUM =
        haskey(config, :REDUNDANT_STEPS_NUM) ? config[:REDUNDANT_STEPS_NUM] : 100
    max_sim_time = haskey(config, :max_sim_time) ? config[:max_sim_time] : Inf
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
        ADAPT_COEFFI_PS,
        ADAPT_COEFFI_VS_GLOBAL,
        ADAPT_COEFFI_VS_LOCAL,
        ADAPT_VS_MODE,
        ADAPT_COEFFI_VS_LOHNER,
        TOLERANCE,
        ST_CHECK_INTERVAL,
        REDUNDANT_STEPS_NUM,
        max_sim_time,
    )
end
function Solver(; kwargs...)
    CFL = haskey(kwargs, :CFL) ? kwargs[:CFL] : 0.4
    AMR_DYNAMIC_PS_MAXLEVEL =
        haskey(kwargs, :AMR_DYNAMIC_PS_MAXLEVEL) ? kwargs[:AMR_DYNAMIC_PS_MAXLEVEL] :
        kwargs[:AMR_PS_MAXLEVEL]
    PS_DYNAMIC_AMR = haskey(kwargs, :PS_DYNAMIC_AMR) ? kwargs[:PS_DYNAMIC_AMR] : true
    VS_DYNAMIC_AMR = haskey(kwargs, :VS_DYNAMIC_AMR) ? kwargs[:VS_DYNAMIC_AMR] : true
    ADAPT_COEFFI_PS = haskey(kwargs, :ADAPT_COEFFI_PS) ? kwargs[:ADAPT_COEFFI_PS] : 0.25
    ADAPT_COEFFI_VS_GLOBAL =
        haskey(kwargs, :ADAPT_COEFFI_VS_GLOBAL) ? kwargs[:ADAPT_COEFFI_VS_GLOBAL] : 0.125
    ADAPT_COEFFI_VS_LOCAL =
        haskey(kwargs, :ADAPT_COEFFI_VS_LOCAL) ? kwargs[:ADAPT_COEFFI_VS_LOCAL] : 1e-2
    ADAPT_VS_MODE =
        haskey(kwargs, :ADAPT_VS_MODE) ? Symbol(kwargs[:ADAPT_VS_MODE]) : :lohner
    ADAPT_COEFFI_VS_LOHNER =
        haskey(kwargs, :ADAPT_COEFFI_VS_LOHNER) ? kwargs[:ADAPT_COEFFI_VS_LOHNER] : 0.6
    TOLERANCE = haskey(kwargs, :TOLERANCE) ? kwargs[:TOLERANCE] : 1e-6
    ST_CHECK_INTERVAL =
        haskey(kwargs, :ST_CHECK_INTERVAL) ? kwargs[:ST_CHECK_INTERVAL] : 100
    REDUNDANT_STEPS_NUM =
        haskey(kwargs, :REDUNDANT_STEPS_NUM) ? kwargs[:REDUNDANT_STEPS_NUM] : 100
    max_sim_time = haskey(kwargs, :max_sim_time) ? kwargs[:max_sim_time] : Inf
    return Solver{kwargs[:DIM],kwargs[:NDF]}(
        CFL,
        kwargs[:AMR_PS_MAXLEVEL],
        AMR_DYNAMIC_PS_MAXLEVEL,
        kwargs[:AMR_VS_MAXLEVEL],
        kwargs[:flux],
        kwargs[:time_marching],
        PS_DYNAMIC_AMR,
        VS_DYNAMIC_AMR,
        ADAPT_COEFFI_PS,
        ADAPT_COEFFI_VS_GLOBAL,
        ADAPT_COEFFI_VS_LOCAL,
        ADAPT_VS_MODE,
        ADAPT_COEFFI_VS_LOHNER,
        TOLERANCE,
        ST_CHECK_INTERVAL,
        REDUNDANT_STEPS_NUM,
        max_sim_time,
    )
end
function Solver(solver::Solver{DIM,NDF}; kwargs...) where {DIM,NDF}
    new_args = [get(kwargs, f, getfield(solver, f)) for f in fieldnames(Solver)]
    return Solver{DIM,NDF}(new_args...)
end

"""
$(TYPEDEF)

Container for the optional user-supplied callbacks that steer adaptive mesh refinement.
Construct with `UDF(; static_ps_refine_flag = …, dynamic_ps_refine_flag = …)`; any field left
unset defaults to a no-op. See the [User-defined functions](@ref) page for the full convention
and worked examples.

## Fields

$(TYPEDFIELDS)
"""
mutable struct UDF
    """
    Static, geometry-driven physical-space refinement flag. Evaluated once during the initial
    refinement at setup ([`pre_refine!`](@ref), called inside [`initialize`](@ref)). Signature

        static_ps_refine_flag(midpoint::Vector{Float64}, ds::Vector{Float64},
                              kinfo::KInfo, level::Int) -> Bool

    `midpoint`/`ds` are the candidate cell's centre/size, `level` its current refinement level.
    Return `true` to force the cell to refine; such cells are also protected from coarsening.
    Default (`null_udf`): never force refinement.
    """
    static_ps_refine_flag::Function
    """
    Dynamic, solution-driven physical-space refinement flag. Evaluated every refinement step
    ([`adaptive_mesh_refinement!`](@ref)) to *gate* the Löhner sensor. Signature

        dynamic_ps_refine_flag(ps_data::AbstractPsData, level::Int, ka::KA) -> Bool

    Return `false` to forbid dynamic refinement of `ps_data` (e.g. to freeze the mesh outside a
    region of interest); `true` lets the sensor decide. Default (unset): always allow.
    """
    dynamic_ps_refine_flag::Function
    """
    Reserved for a static velocity-space refinement flag. **Currently unused** — the field is stored but never invoked.
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
    VTK cell type in physical space. Available options: `Pixel` and `Triangle` for 2D; `Voxel` and `Tetra` for 3D. Default is [`Triangle`](@ref) for 2D and [`Tetra`](@ref) for 3D. A `Vector` of cell types (e.g. `[Triangle, Pixel]`) may also be given; in that case a separate output (pvtu for the flow field, and its own animation collection) is written for every cell type, with the cell-type name appended to the corresponding file name.
    """
    vtk_celltype::Union{Type{<:AbstractVTKCellType},Vector}
    """
    VTK cell type in velocity space. Available options: `Pixel` for 2D.
    """
    vs_vtk_celltype::Type{Tv} where {Tv<:AbstractVTKCellType}
    """
    Time interval between animation frames. `0` (default) disables animation output, making [`check_for_animsave!`](@ref) a no-op.
    """
    anim_dt::Float64
    """
    Optional callback selecting which physical cells write their (per-cell) velocity-space
    distribution. Two call conventions, one per output path:

      - Final result ([`save_result`](@ref)): `vs_output_criterion(ps_data, ka) -> Bool` —
        return `true` to include this cell's velocity space. Default (unset): include every cell.
      - Animation frames ([`check_for_animsave!`](@ref)): `vs_output_criterion(; ps_data, ka) -> (id::Int, flag::Bool)` — `flag` selects the cell and `id` names its output file.
        Default (unset): write no per-cell velocity space.
        See the [User-defined functions](@ref) page.
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
            setfield!(output, i, config[i])
        end
    end
    return output
end
function Output(solver::Solver{DIM,NDF}; kwargs...) where {DIM,NDF}
    output = Output(DIM==2 ? Triangle : Tetra, DIM==2 ? Pixel : Voxel, 0.0, null_udf, 0)
    for i in fieldnames(Output)
        if haskey(kwargs, i)
            setfield!(output, i, kwargs[i])
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
struct Configure{DIM,NDF}<:AbstractConfig{DIM,NDF}
    """
    Range of the simulated domain. As an example, for 2D case, it should be aligned as [xmin,xmax,ymin,ymax].
    """
    geometry::Vector{Float64}
    """
    Number of tree roots for each dimension. For 2D case, it should be aligned as [x_num,y_num]
    """
    trees_num::Vector{Int64}
    """
    For `Vector` type, it represents the range of the velocity space. For 2D case, it should be
    aligned as `[umin,umax,vmin,vmax]`. **Implicit requirement: the origin `v = 0` must lie on a
    velocity-grid corner**, i.e. for every dimension `(0 - min)/(max - min)*vs_trees_num` must be
    an integer (checked at construction by [`check_vs_setting`](@ref)). This keeps the upwinding
    split at `v·n = 0` unambiguous and is easy to satisfy — e.g. a symmetric range `[-a, a]` with
    an even `vs_trees_num`. The `AbstractQuadrature` type is reserved for other quadrature rules
    like Gauss-Hermite.
    """
    quadrature::Union{AbstractQuadrature,Vector{Float64}}
    """
    Number of tree roots for each dimension in velocity space. Together with `quadrature` it must place the origin `v = 0` on a root-grid corner (see `quadrature`; checked by [`check_vs_setting`](@ref)).
    """
    vs_trees_num::Vector{Int64}
    """
    Initial condition defined by [`AbstractInitCond`](@ref).
    """
    IC::AbstractInitCond
    """
    Types of the domain boundary. For 2D case, the vector should catain 4 elements corresponding to the 4 domain boundaries. The element type is defined by [`Domain`](@ref).
    """
    domain::Vector{Domain}
    """
    Immersed boundary. Multiple boundaries are supported. The element type is defined by [`AbstractBoundary`](@ref)
    """
    IB::Vector{AbstractBoundary}
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
"""
$(TYPEDSIGNATURES)
Validate the Cartesian velocity-space setting (`quadrature` range + `vs_trees_num`).

The kinetic scheme carries an **implicit requirement: the origin `v = 0` must lie on a
velocity-grid corner** (a root-cell vertex), never inside a cell. This keeps the upwinding split
at `v·n = 0` unambiguous (no cell whose centre is exactly `v = 0`) and, because refinement only
bisects cells, guarantees 0 stays on a corner at every level. Equivalently, for each dimension
`(0 - min)/(max - min) * vs_trees_num` must be an integer in `0:vs_trees_num`.

The requirement is left to the user to satisfy through a sensible choice of range and root count —
it is trivial, e.g. a symmetric range `[-a, a]` with an even `vs_trees_num`. This function only
checks it, throwing an `ErrorException` that names the offending dimension when violated. It is a
no-op for a non-Cartesian `AbstractQuadrature`.
"""
function check_vs_setting(quadrature, vs_trees_num, DIM::Integer)
    quadrature isa AbstractVector || return nothing   # non-Cartesian quadrature: requirement N/A
    length(quadrature) == 2DIM || error(
        "`quadrature` must have length 2*DIM = $(2DIM) (got $(length(quadrature))): an increasing [min,max] pair per velocity dimension.",
    )
    length(vs_trees_num) == DIM ||
        error("`vs_trees_num` must have length DIM = $DIM (got $(length(vs_trees_num))).")
    for d = 1:DIM
        lo = quadrature[2d-1];
        hi = quadrature[2d];
        n = vs_trees_num[d]
        hi > lo || error(
            "velocity-space dimension $d: `quadrature` range [$lo, $hi] must be increasing.",
        )
        n >= 1 || error("velocity-space dimension $d: `vs_trees_num[$d]` = $n must be ≥ 1.")
        k = -lo * n / (hi - lo)        # fractional root-grid index of the v = 0 vertex
        i = round(k)
        (0 <= i <= n && abs(k - i) <= 1e-8 * max(1, n)) || error(
            "Invalid velocity-space setting in dimension $d: the origin v = 0 must lie on a " *
            "root-grid corner. With range [$lo, $hi] and `vs_trees_num[$d]` = $n, 0 falls at " *
            "fractional cell index $k (it must be an integer in 0:$n). Adjust the `quadrature` " *
            "range and/or `vs_trees_num` so that (0 - min)/(max - min)*vs_trees_num is an integer " *
            "— e.g. a symmetric range [-a, a] with an even `vs_trees_num`.",
        )
    end
    return nothing
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
    check_vs_setting(config[:quadrature], config[:vs_trees_num], config[:DIM])
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
    check_vs_setting(kwargs[:quadrature], kwargs[:vs_trees_num], DIM)
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
function Configure(solver::Solver, config::Configure{DIM,NDF}; kwargs...) where {DIM,NDF}
    return Configure(
        solver;
        geometry = haskey(kwargs, :geometry) ? kwargs[:geometry] : config.geometry,
        trees_num = haskey(kwargs, :trees_num) ? kwargs[:trees_num] : config.trees_num,
        quadrature = haskey(kwargs, :quadrature) ? kwargs[:quadrature] : config.quadrature,
        vs_trees_num = haskey(kwargs, :vs_trees_num) ? kwargs[:vs_trees_num] :
                       config.vs_trees_num,
        IC = haskey(kwargs, :IC) ? kwargs[:IC] : config.IC,
        domain = haskey(kwargs, :domain) ? kwargs[:domain] : config.domain,
        IB = haskey(kwargs, :IB) ? kwargs[:IB] : config.IB,
        gas = haskey(kwargs, :gas) ? kwargs[:gas] : config.gas,
        output = haskey(kwargs, :output) ? kwargs[:output] : config.output,
        user_defined = haskey(kwargs, :user_defined) ? kwargs[:user_defined] :
                       config.user_defined,
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
    return Residual(ones(DIM+2), zeros(DIM+2), zeros(DIM+2), 0)
end

"""
$(TYPEDEF)

Structure of simulation status. It contains the real time information, and is updated every step.

## Fields

$(TYPEDFIELDS)
"""
mutable struct Status
    """
    Maximum absolute value of the gradients of conserved variables.
    """
    gradmax::Vector{Float64}
    """
    Maximum number of velocity cells in a physical cell.
    """
    max_vs_num::Int
    """
    Total number of phase grids.
    """
    total_phase_num::Int
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
    Number of steps that have been marched.
    """
    step::Int
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
    """
    Reusable vector of MPI requests for asynchronous communication.
    """
    mpi_reqs::Vector{MPI.Request}
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
        zeros(DIM+2),
        0,
        0,
        Δt_ξ,
        Δt_ξ,
        0.0,
        0,
        1,
        1,
        1,
        Residual(DIM),
        Ref(false),
        MPI.Request[],
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
        zeros(DIM+2),
        0,
        0,
        Δt_ξ,
        Δt_ξ,
        0.0,
        0,
        1,
        1,
        1,
        Residual(DIM),
        Ref(false),
        MPI.Request[],
    )
end

"""
$(TYPEDEF)
Structure of required information during a simulation.
$(TYPEDFIELDS)
"""
mutable struct KInfo{DIM,NDF}
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
end
function KInfo(config::Dict)
    return KInfo{config[:DIM],config[:NDF]}(
        Configure(config),
        Forest(config[:DIM]),
        Status(config),
    )
end
function KInfo(config::Configure{DIM,NDF}) where {DIM,NDF}
    return KInfo{DIM,NDF}(config, Forest(DIM), Status(config))
end
"""
$(TYPEDSIGNATURES)
Update information in [`KInfo`](@ref) to maintain self consistency.
"""
function KInfo!(kinfo::KInfo{DIM,NDF}; kwargs...) where {DIM,NDF}
    if haskey(kwargs, :solver)
        solver = kwargs[:solver]
        config = Configure(solver, kinfo.config)
    else
        config = haskey(kwargs, :config) ? kwargs[:config] : kinfo.config
    end
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
    kinfo.status.Δt = Δt_ξ;
    kinfo.status.Δt_ξ = Δt_ξ
    kinfo.status.residual.redundant_step = 0
    kinfo.status.residual.residual = Inf .* ones(DIM+2)
    kinfo.config = config
end

mutable struct P4estPsData
    ps_data::Ptr{Nothing}
end

"""
$(TYPEDEF)
Julia-managed MPI buffers for the ghost layer.  All memory is owned by the Julia
GC; no manual `sc_free` is required.
$(TYPEDFIELDS)
"""
mutable struct GhostBuffer
    """
    Flat receive buffer for ghost cell data (df, w, …).
    """
    ghost_datas::Vector{Float64}
    """
    Flat receive buffer for ghost cell slopes (sdf, sw).
    """
    ghost_slopes::Vector{Float64}
    """
    Flat receive buffer for ghost VS structure (weight, level, midpoint).
    """
    ghost_structures::Vector{Float64}
    """
    Per-mirror send buffer for data, one entry per mirror quadrant.
    """
    mirror_data_bufs::Vector{Vector{Float64}}
    """
    Per-mirror send buffer for slopes, one entry per mirror quadrant.
    """
    mirror_slope_bufs::Vector{Vector{Float64}}
    """
    Per-mirror send buffer for VS structure, one entry per mirror quadrant.
    """
    mirror_structure_bufs::Vector{Vector{Float64}}
end

"""
$(TYPEDEF)
Julia-managed metadata for the ghost layer: vs_nums, flat-buffer offsets, and
pre-computed element-count arrays that are reused every time-step without reallocation.
$(TYPEDFIELDS)
"""
struct GhostInfo
    """
    vs_num of each ghost quadrant (length = n_ghosts).
    """
    ghost_vsnums::Vector{Int}
    """
    vs_num of each mirror quadrant (length = n_mirrors).
    """
    mirror_vsnums::Vector{Int}
    """
    Float64-element offset into ghost_datas for each ghost quadrant.
    """
    ghost_data_offsets::Vector{Int}
    """
    Float64-element offset into ghost_slopes for each ghost quadrant.
    """
    ghost_slope_offsets::Vector{Int}
    """
    Float64-element offset into ghost_structures for each ghost quadrant.
    """
    ghost_structure_offsets::Vector{Int}
    """
    Cached Float64-element count of ghost data buffer per ghost quadrant.
    """
    ghost_data_szs::Vector{Int}
    """
    Cached Float64-element count of ghost slope buffer per ghost quadrant.
    """
    ghost_slope_szs::Vector{Int}
    """
    Cached Float64-element count of ghost data buffer per mirror quadrant.
    """
    mirror_data_szs::Vector{Int}
    """
    Cached Float64-element count of ghost slope buffer per mirror quadrant.
    """
    mirror_slope_szs::Vector{Int}
    """
    Physical-space refinement level of each ghost quadrant.
    """
    ghost_levels::Vector{Int8}
    """
    Physical-space refinement level of each mirror quadrant.
    """
    mirror_levels::Vector{Int8}
end

"""
$(TYPEDEF)
Structure of cells data managed by `Julia`.
$(TYPEDFIELDS)
"""
mutable struct PsTrees{DIM,NDF}
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
    ghost_buffer::GhostBuffer
    ghost_wrap::Vector{AbstractGhostPsData}
    ghost_info::GhostInfo
end

"""
$(TYPEDEF)
Structure of fields data.
$(TYPEDFIELDS)
"""
mutable struct Field{DIM,NDF}
    """
    Cells data defined by [`PsTrees`](@ref).
    """
    trees::PsTrees{DIM,NDF}
    """
    Mapping between faces and cells. The element type is defined by [`AbstractFace`](@ref).
    """
    faces::Vector{AbstractFace}
    """
    Pre-collected immersed-boundary data (donor cells, solid cells, IB faces).
    """
    immersed_boundaries::Vector{ImmersedBoundary{DIM,NDF}}
end

"""
$(TYPEDEF)
Structure of computed data.
$(TYPEDFIELDS)
"""
mutable struct KData{DIM,NDF}
    field::Field{DIM,NDF}
    ghost::Ghost
    KData(trees::PsTrees{DIM,NDF}) where {DIM,NDF} = (
        n = new{DIM,NDF}();
        n.field = Field{DIM,NDF}(
            trees,
            Vector{AbstractFace}(undef, 0),
            ImmersedBoundary{DIM,NDF}[],
        );
        n
    )
end

"""
$(TYPEDEF)
Structure that collects all of the data.
$(TYPEDFIELDS)
"""
mutable struct KA{DIM,NDF}
    kinfo::KInfo{DIM,NDF}
    kdata::KData{DIM,NDF}
end
# NB: with the fields typed `KInfo{DIM,NDF}`/`KData{DIM,NDF}`, Julia auto-generates
# the outer constructor `KA(::KInfo{DIM,NDF}, ::KData{DIM,NDF})` that infers {DIM,NDF};
# an explicit one would be a duplicate (a hard precompile error on Julia ≥1.12).
# partition
struct TransferData{DIM,NDF}
    encs::Vector{Int}
    w::Vector{Float64}
    vs_levels::Vector{Int8}
    vs_midpoints::Vector{Float64}
    vs_df::Vector{Float64}
end
function TransferData(DIM::Integer, NDF::Integer, ps_num::Integer, total_vs_num::Integer)
    return TransferData{DIM,NDF}(
        Vector{Int}(undef, (SOLID_CELL_ID_NUM+1)*ps_num),
        Vector{Float64}(undef, (DIM+2) * ps_num),
        Vector{Int8}(undef, total_vs_num),
        Vector{Float64}(undef, DIM * total_vs_num),
        Vector{Float64}(undef, NDF * total_vs_num),
    )
end

mutable struct TransferInit
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
