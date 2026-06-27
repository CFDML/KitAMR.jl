module R2DVJEffi

using KitAMR
using MPI
using JLD2
using Printf

const MODES = ("both_on", "both_off", "both_on_3d", "both_off_3d")

const OPTS = Dict{String,String}()
for (k, v) in ENV
    startswith(k, "VJE_") && (OPTS[k] = v)
end
for a in ARGS
    occursin('=', a) || continue
    k, v = split(a, '=', limit = 2)
    startswith(k, "VJE_") && (OPTS[k] = v)
end

envs(k, d) = get(OPTS, k, d)
envi(k, d) = parse(Int, get(OPTS, k, string(d)))
envf(k, d) = parse(Float64, get(OPTS, k, string(d)))
envb(k, d) = lowercase(get(OPTS, k, d ? "1" : "0")) in ("1", "true", "yes", "on")

function env_interval(k, d)
    s = lowercase(envs(k, d))
    s == "auto" && return :auto
    value = parse(Int, s)
    value > 0 || error("$k must be a positive integer or auto, got $s")
    return value
end

function with_mpi(f::Function)
    started_here = !MPI.Initialized()
    started_here && MPI.Init()
    try
        return f()
    finally
        started_here && !MPI.Finalized() && MPI.Finalize()
    end
end

const DIM = envi("VJE_DIM", 2)
DIM in (2, 3) || error("VJE_DIM must be 2 or 3, got $DIM")
const NDF = DIM == 2 ? 2 : 1
const LOG_TAG = DIM == 2 ? "r2d_vj_effi" : "r3d_vj_effi"

const PS_BASE = envi("VJE_PS_BASE", 8)
const PS_BASE_X = haskey(OPTS, "VJE_PS_BASE_X") ? envi("VJE_PS_BASE_X", PS_BASE) : PS_BASE
const PS_BASE_Y = haskey(OPTS, "VJE_PS_BASE_Y") ? envi("VJE_PS_BASE_Y", 2 * PS_BASE_X) : 2 * PS_BASE_X
const PS_BASE_Z = haskey(OPTS, "VJE_PS_BASE_Z") ? envi("VJE_PS_BASE_Z", PS_BASE_Y) : PS_BASE_Y
const PS_LMAX = envi("VJE_PS_LMAX", 5)
const VS_BASE = envi("VJE_VS_BASE", 8)
const CFL = envf("VJE_CFL", 0.4)
const MAX_T = envf("VJE_MAX_SIM_TIME", 0.1)
const ANIM_DT = envf("VJE_ANIM_DT", 0.025)
const MAX_STEP = haskey(OPTS, "VJE_MAX_STEPS") ? envi("VJE_MAX_STEPS", 0) : typemax(Int)
const OUTDIR = envs("VJE_OUTDIR", joinpath(@__DIR__, "out"))
const PS_INTERVAL = env_interval("VJE_PS_INTERVAL", "auto")
const VS_INTERVAL = env_interval("VJE_VS_INTERVAL", "auto")
const PARTITION_INTERVAL = env_interval("VJE_PARTITION_INTERVAL", "auto")
const AUTO_PS_TRAVEL_FRACTION = envf("VJE_AUTO_PS_TRAVEL_FRACTION", 2.0)

const RHO_JET = envf("VJE_RHO_JET", 1.0)
const U_JET = envf("VJE_U_JET", 6.0)
const V_JET = envf("VJE_V_JET", 0.0)
const W_JET = envf("VJE_W_JET", 0.0)
const LAMBDA_JET = envf("VJE_LAMBDA_JET", 0.5)
const RHO_VAC = envf("VJE_RHO_VAC", 1.0e-4)
const LAMBDA_VAC = envf("VJE_LAMBDA_VAC", 2.0)
const KN = envf("VJE_KN", 1.0e-3)
const JET_RADIUS = envf("VJE_JET_RADIUS", 0.16)
const JET_EDGE_W = envf("VJE_JET_EDGE_WIDTH", 0.015)

const GAMMA = 5 / 3
const THETA_JET = 1.0 / (2.0 * LAMBDA_JET)
const THETA_VAC = 1.0 / (2.0 * LAMBDA_VAC)
const THETA_MIN = min(THETA_JET, THETA_VAC)
const THETA_STAT = max(THETA_JET, THETA_VAC)
const SIGMA_MIN = sqrt(THETA_MIN)
const U_MAX = DIM == 2 ? hypot(U_JET, V_JET) : sqrt(U_JET^2 + V_JET^2 + W_JET^2)
const THETA_MAX = THETA_STAT + 0.5 * U_MAX^2 * (GAMMA - 1.0) / GAMMA
const SIGMA_MAX = sqrt(THETA_MAX)
const VRANGE = haskey(OPTS, "VJE_VRANGE") ? envf("VJE_VRANGE", 0.0) : U_MAX + 5.0 * SIGMA_MAX
const MIN_VS_PER_DIM = ceil(Int, 6.0 * VRANGE / SIGMA_MIN)
const VS_LMAX = haskey(OPTS, "VJE_VS_LMAX") ? envi("VJE_VS_LMAX", 0) :
                max(0, ceil(Int, log2(MIN_VS_PER_DIM / VS_BASE)))

# Computational box: the jet enters from the left-boundary center. Keep every non-inlet
# boundary at least one fastest-molecule travel distance away from the inlet point:
# downstream length L and transverse span 2L. y/z roots default to twice x-roots to keep
# physical-space cells cubic at the root and finest levels.
const DOMAIN_TRAVEL_TIME = envf("VJE_DOMAIN_TRAVEL_TIME", 0.1)
const DOMAIN_L = VRANGE * DOMAIN_TRAVEL_TIME
const GEOM = DIM == 2 ? [0.0, DOMAIN_L, -DOMAIN_L, DOMAIN_L] :
                         [0.0, DOMAIN_L, -DOMAIN_L, DOMAIN_L, -DOMAIN_L, DOMAIN_L]
const XC0 = GEOM[1]
const XC1 = GEOM[2]
const YC0 = GEOM[3]
const YC1 = GEOM[4]
const YMID = 0.5 * (YC0 + YC1)
const ZC0 = DIM == 3 ? GEOM[5] : 0.0
const ZC1 = DIM == 3 ? GEOM[6] : 0.0
const ZMID = 0.5 * (ZC0 + ZC1)

const REF_PS_X = PS_BASE_X * 2^PS_LMAX
const REF_PS_Y = PS_BASE_Y * 2^PS_LMAX
const REF_PS_Z = PS_BASE_Z * 2^PS_LMAX
const REF_VS_PER_DIM = VS_BASE * 2^VS_LMAX
const REF_PS_CELLS = DIM == 2 ? REF_PS_X * REF_PS_Y : REF_PS_X * REF_PS_Y * REF_PS_Z
const REF_VS_CELLS = REF_VS_PER_DIM^DIM
const DENSE_PHASE_EST = Float64(REF_PS_CELLS) * Float64(REF_VS_CELLS)
const ALLOW_DENSE = envb("VJE_ALLOW_DENSE", false)
const MAX_UNIFORM_PHASE = envf("VJE_MAX_UNIFORM_PHASE", 3.0e8)
const BOTH_OFF_ROOTS_PER_RANK = envf("VJE_BOTH_OFF_ROOTS_PER_RANK", 1.0)
const BOTH_OFF_MIN_ROOTS = envi("VJE_BOTH_OFF_MIN_ROOTS", 0)

function check_inputs()
    PS_BASE > 0 || error("VJE_PS_BASE must be positive")
    PS_BASE_X > 0 || error("VJE_PS_BASE_X must be positive")
    PS_BASE_Y > 0 || error("VJE_PS_BASE_Y must be positive")
    PS_BASE_Z > 0 || error("VJE_PS_BASE_Z must be positive")
    PS_LMAX >= 0 || error("VJE_PS_LMAX must be non-negative")
    VS_BASE > 0 || error("VJE_VS_BASE must be positive")
    VS_LMAX >= 0 || error("VJE_VS_LMAX must be non-negative")
    BOTH_OFF_ROOTS_PER_RANK > 0 || error("VJE_BOTH_OFF_ROOTS_PER_RANK must be positive")
    BOTH_OFF_MIN_ROOTS >= 0 || error("VJE_BOTH_OFF_MIN_ROOTS must be non-negative")
    RHO_JET > 0 || error("VJE_RHO_JET must be positive")
    RHO_VAC > 0 || error("VJE_RHO_VAC must be positive")
    LAMBDA_JET > 0 || error("VJE_LAMBDA_JET must be positive")
    LAMBDA_VAC > 0 || error("VJE_LAMBDA_VAC must be positive")
    JET_RADIUS > 0 || error("VJE_JET_RADIUS must be positive")
    JET_EDGE_W > 0 || error("VJE_JET_EDGE_WIDTH must be positive")
    AUTO_PS_TRAVEL_FRACTION > 0 || error("VJE_AUTO_PS_TRAVEL_FRACTION must be positive")
    DOMAIN_TRAVEL_TIME > 0 || error("VJE_DOMAIN_TRAVEL_TIME must be positive")
    (YMID - JET_RADIUS > YC0 && YMID + JET_RADIUS < YC1) ||
        error("jet radius must fit inside y-domain")
    if DIM == 3
        (ZMID - JET_RADIUS > ZC0 && ZMID + JET_RADIUS < ZC1) ||
            error("jet radius must fit inside z-domain")
    end
    dx_root = (XC1 - XC0) / PS_BASE_X
    dy_root = (YC1 - YC0) / PS_BASE_Y
    isapprox(dx_root, dy_root; rtol = 1.0e-10, atol = 1.0e-12) ||
        error("physical root cells must be square; got dx=$dx_root dy=$dy_root. " *
              "Use VJE_PS_BASE_Y = 2 * VJE_PS_BASE_X for the default box.")
    if DIM == 3
        dz_root = (ZC1 - ZC0) / PS_BASE_Z
        isapprox(dx_root, dz_root; rtol = 1.0e-10, atol = 1.0e-12) ||
            error("physical root cells must be cubic; got dx=$dx_root dz=$dz_root. " *
                  "Use VJE_PS_BASE_Z = 2 * VJE_PS_BASE_X for the default box.")
    end
    return nothing
end

check_inputs()

function smooth_aperture_profile(s, center, radius, edge_width)
    return clamp(0.5 * (1.0 - tanh((abs(s - center) - radius) / edge_width)), 0.0, 1.0)
end

jet_profile_y(y) = smooth_aperture_profile(y, YMID, JET_RADIUS, JET_EDGE_W)
jet_profile_yz(y, z) = clamp(0.5 * (1.0 - tanh((hypot(y - YMID, z - ZMID) -
    JET_RADIUS) / JET_EDGE_W)), 0.0, 1.0)

function jet_primitive(a)
    rho = RHO_VAC + a * (RHO_JET - RHO_VAC)
    ux = a * U_JET
    uy = a * V_JET
    theta = THETA_VAC + a * (THETA_JET - THETA_VAC)
    lambda = 1.0 / (2.0 * theta)
    if DIM == 2
        return [rho, ux, uy, lambda]
    else
        uz = a * W_JET
        return [rho, ux, uy, uz, lambda]
    end
end

left_jet_primitive(a) = jet_primitive(a)
vacuum_jet_init(midpoint, kinfo) = left_jet_primitive(0.0)
left_aperture_weight(midpoint) = DIM == 2 ? jet_profile_y(midpoint[2]) :
    jet_profile_yz(midpoint[2], midpoint[3])
left_composite_weights(midpoint) = (a = left_aperture_weight(midpoint); (a, 1.0 - a))
const LEFT_JET_BC_PRIM = left_jet_primitive(1.0)
const LEFT_WALL_BC_PRIM = DIM == 2 ? [RHO_VAC, 0.0, 0.0, LAMBDA_VAC] :
                                     [RHO_VAC, 0.0, 0.0, 0.0, LAMBDA_VAC]

no_vs(ps_data, ka) = false
no_vs(; ps_data, ka) = (0, false)

function temperature_lohner_criterion(ps_data, level, ka)
    value = 0.0
    temperature_row = size(ps_data.lohner, 1)
    @inbounds for dir in axes(ps_data.lohner, 2)
        value = max(value, ps_data.lohner[temperature_row, dir])
    end
    return value
end

function static_reference_ps_refine(midpoint, ds, kinfo, level)
    return !kinfo.config.solver.AMR_PS_DYNAMIC &&
           level < kinfo.config.solver.AMR_PS_MAXLEVEL
end

mpi_size_or_one() =
    MPI.Initialized() && !MPI.Finalized() ? MPI.Comm_size(MPI.COMM_WORLD) : 1

function both_off_min_root_trees()
    by_rank = ceil(Int, BOTH_OFF_ROOTS_PER_RANK * mpi_size_or_one())
    return max(BOTH_OFF_MIN_ROOTS, by_rank, 1)
end

function static_reference_layout(ref_dims::NTuple{ND,Int}, maxlevel::Integer,
                                 min_root_trees::Integer) where {ND}
    min_root_trees <= prod(ref_dims) ||
        error("both_off needs at least $min_root_trees root trees, but the finest " *
              "physical mesh only has $(prod(ref_dims)) cells.")

    static_level = 0
    for level in 0:maxlevel
        divisor = 2^level
        all(d -> ref_dims[d] % divisor == 0, 1:ND) || break
        root_dims = ntuple(d -> ref_dims[d] ÷ divisor, Val(ND))
        prod(root_dims) >= min_root_trees || break
        static_level = level
    end
    root_dims = ntuple(d -> ref_dims[d] ÷ 2^static_level, Val(ND))
    return root_dims, static_level
end

struct ModeSpec
    mode::String
    dimension::Int
    ps_dynamic::Bool
    vs_dynamic::Bool
    ps_trees_x::Int
    ps_trees_y::Int
    ps_trees_z::Int
    vs_trees::Int
    ps_amr_level::Int
    vs_amr_level::Int
    prerefine_steps::Int
end

function mode_spec(mode::AbstractString)
    mode in MODES || error("mode must be one of $(join(MODES, ", ")), got $mode")
    mode_is_3d = endswith(mode, "_3d")
    mode_is_3d && DIM == 2 && error("mode=$mode requires VJE_DIM=3")
    !mode_is_3d && DIM == 3 && (mode = string(mode, "_3d"))
    base_mode = replace(String(mode), "_3d" => "")
    ps_dynamic = base_mode == "both_on"
    vs_dynamic = base_mode == "both_on"
    if ps_dynamic
        ps_trees_x = PS_BASE_X
        ps_trees_y = PS_BASE_Y
        ps_trees_z = PS_BASE_Z
        ps_amr = PS_LMAX
        prerefine = PS_LMAX
    elseif DIM == 2
        root_dims, ps_amr = static_reference_layout(
            (REF_PS_X, REF_PS_Y), PS_LMAX, both_off_min_root_trees())
        ps_trees_x, ps_trees_y = root_dims
        ps_trees_z = 1
        prerefine = 0
    else
        root_dims, ps_amr = static_reference_layout(
            (REF_PS_X, REF_PS_Y, REF_PS_Z), PS_LMAX, both_off_min_root_trees())
        ps_trees_x, ps_trees_y, ps_trees_z = root_dims
        prerefine = 0
    end
    vs_trees = vs_dynamic ? VS_BASE : VS_BASE * 2^VS_LMAX
    vs_amr = vs_dynamic ? VS_LMAX : 0
    if base_mode != "both_on" && DENSE_PHASE_EST > MAX_UNIFORM_PHASE && !ALLOW_DENSE
        error("mode=$mode estimates $(DENSE_PHASE_EST) dense phase cells; lower VJE_PS_LMAX/VJE_VS_LMAX or set VJE_ALLOW_DENSE=1 after checking memory")
    end
    return ModeSpec(String(mode), DIM, ps_dynamic, vs_dynamic, ps_trees_x, ps_trees_y,
                    ps_trees_z, vs_trees, ps_amr, vs_amr, prerefine)
end

function build_config(spec::ModeSpec)
    solver = Solver(;
        DIM = DIM,
        NDF = NDF,
        CFL = CFL,
        AMR_PS_MAXLEVEL = spec.ps_amr_level,
        AMR_VS_MAXLEVEL = spec.vs_amr_level,
        AMR_PS_DYNAMIC = spec.ps_dynamic,
        AMR_VS_DYNAMIC = spec.vs_dynamic,
        AUTO_AMR_PS_TRAVEL_FRACTION = AUTO_PS_TRAVEL_FRACTION,
        flux = CAIDVM,
        time_marching = CIP_Marching,
        max_sim_time = MAX_T,
    )
    udf = UDF(;
        static_ps_refine_flag = static_reference_ps_refine,
        dynamic_ps_adapt_criterion = temperature_lohner_criterion,
    )

    gas = Gas(; K = DIM == 2 ? 1.0 : 0.0, Kn = KN, ω = 0.81, ωᵣ = 0.81)
    output = Output(solver;
        vtk_celltype = DIM == 2 ? [Pixel, Triangle] : [Voxel, Tetra],
        vs_vtk_celltype = DIM == 2 ? Pixel : Voxel,
        anim_dt = ANIM_DT,
        vs_output_criterion = no_vs,
    )
    left_domain = Domain(Composite, 1, CompositeBC(left_composite_weights,
        Domain(SuperSonicInflow, 1, LEFT_JET_BC_PRIM),
        Domain(Maxwellian, 1, LEFT_WALL_BC_PRIM)))

    if DIM == 2
        return Configure(solver;
            geometry = GEOM,
            trees_num = [spec.ps_trees_x, spec.ps_trees_y],
            quadrature = [-VRANGE, VRANGE, -VRANGE, VRANGE],
            vs_trees_num = [spec.vs_trees, spec.vs_trees],
            IC = PCoordFn(vacuum_jet_init),
            domain = [left_domain, Domain(UniformOutflow, 2),
                      Domain(UniformOutflow, 3), Domain(UniformOutflow, 4)],
            output = output,
            gas = gas,
            user_defined = udf,
        )
    else
        return Configure(solver;
            geometry = GEOM,
            trees_num = [spec.ps_trees_x, spec.ps_trees_y, spec.ps_trees_z],
            quadrature = [-VRANGE, VRANGE, -VRANGE, VRANGE, -VRANGE, VRANGE],
            vs_trees_num = [spec.vs_trees, spec.vs_trees, spec.vs_trees],
            IC = PCoordFn(vacuum_jet_init),
            domain = [left_domain, Domain(UniformOutflow, 2),
                      Domain(UniformOutflow, 3), Domain(UniformOutflow, 4),
                      Domain(UniformOutflow, 5), Domain(UniformOutflow, 6)],
            output = output,
            gas = gas,
            user_defined = udf,
        )
    end
end

function print_header(spec::ModeSpec)
    rank = MPI.Comm_rank(MPI.COMM_WORLD)
    rank == 0 || return nothing
    root_trees = DIM == 2 ? spec.ps_trees_x * spec.ps_trees_y :
                 spec.ps_trees_x * spec.ps_trees_y * spec.ps_trees_z
    leaves_per_root = 2^(DIM * spec.ps_amr_level)
    roots_per_rank = root_trees / MPI.Comm_size(MPI.COMM_WORLD)
    println("="^78)
    @printf("[%s] mode=%s dim=%d np=%d\n", LOG_TAG, spec.mode, spec.dimension,
            MPI.Comm_size(MPI.COMM_WORLD))
    if DIM == 2
        @printf("[%s] geometry=[%.3f, %.3f] x [%.3f, %.3f]\n",
                LOG_TAG, GEOM[1], GEOM[2], GEOM[3], GEOM[4])
        @printf("[%s] PS trees=(%d,%d) AMR_L=%d dynamic=%s finest=%d x %d\n",
                LOG_TAG, spec.ps_trees_x, spec.ps_trees_y, spec.ps_amr_level,
                string(spec.ps_dynamic), REF_PS_X, REF_PS_Y)
    else
        @printf("[%s] geometry=[%.3f, %.3f] x [%.3f, %.3f] x [%.3f, %.3f]\n",
                LOG_TAG, GEOM[1], GEOM[2], GEOM[3], GEOM[4], GEOM[5], GEOM[6])
        @printf("[%s] PS trees=(%d,%d,%d) AMR_L=%d dynamic=%s finest=%d x %d x %d\n",
                LOG_TAG, spec.ps_trees_x, spec.ps_trees_y, spec.ps_trees_z,
                spec.ps_amr_level, string(spec.ps_dynamic), REF_PS_X, REF_PS_Y, REF_PS_Z)
    end
    @printf("[%s] PS root_trees=%d roots/rank=%.3f leaves/root=%d static_reference=%s\n",
            LOG_TAG, root_trees, roots_per_rank, leaves_per_root,
            string(!spec.ps_dynamic && spec.ps_amr_level > 0))
    @printf("[%s] domain L=%.4g from vrange %.4g over Δt=%.4g\n",
            LOG_TAG,
            DOMAIN_L, VRANGE, DOMAIN_TRAVEL_TIME)
    @printf("[%s] VS trees=%d AMR_L=%d dynamic=%s finest=%d/dim\n",
            LOG_TAG,
            spec.vs_trees, spec.vs_amr_level, string(spec.vs_dynamic), REF_VS_PER_DIM)
    @printf("[%s] max_t=%.4g anim_dt=%.4g vrange=+/-%.4g Kn=%.4g\n",
            LOG_TAG,
            MAX_T, ANIM_DT, VRANGE, KN)
    @printf("[%s] VS estimate: sigma_min=%.4g min_grids/dim=%d (Δv=%.4g, sigma/Δv=%.3f)\n",
            LOG_TAG,
            SIGMA_MIN, MIN_VS_PER_DIM, 2.0 * VRANGE / REF_VS_PER_DIM,
            SIGMA_MIN / (2.0 * VRANGE / REF_VS_PER_DIM))
    @printf("[%s] jet boundary: %s radius=%.4g edge_width=%.4g\n",
            LOG_TAG, DIM == 2 ? "left aperture" : "left circular pipe", JET_RADIUS, JET_EDGE_W)
    @printf("[%s] PS criterion: temperature-only Lohner UDF\n", LOG_TAG)
    @printf("[%s] AMR intervals: PS=%s VS=%s partition=%s vs_balance=true\n",
            LOG_TAG,
            string(PS_INTERVAL), string(VS_INTERVAL), string(PARTITION_INTERVAL))
    @printf("[%s] outdir=%s\n", LOG_TAG, OUTDIR)
    println("="^78)
    return nothing
end

function collect_field(ka)
    nloc = 0
    for tree in ka.kdata.field.trees.data, ps in tree
        isa(ps, KitAMR.InsideSolidData) && continue
        ps.bound_enc < 0 && continue
        nloc += 1
    end

    rows = Matrix{Float64}(undef, nloc, 3 * DIM + 2)
    i = 0
    phase = 0
    maxvs = 0
    mass = 0.0
    energy = 0.0
    minrho = Inf
    maxrho = 0.0
    maxspeed = 0.0

    for tree in ka.kdata.field.trees.data, ps in tree
        isa(ps, KitAMR.InsideSolidData) && continue
        ps.bound_enc < 0 && continue
        i += 1
        volume = prod(ps.ds)
        rho = ps.prim[1]
        theta = 1.0 / (2.0 * ps.prim[end])
        speed = sqrt(sum(ps.prim[d + 1]^2 for d in 1:DIM))
        @inbounds for d in 1:DIM
            rows[i, d] = ps.midpoint[d]
            rows[i, DIM + d] = ps.ds[d]
            rows[i, 2 * DIM + 1 + d] = ps.prim[d + 1]
        end
        rows[i, 2 * DIM + 1] = rho
        rows[i, 3 * DIM + 2] = theta
        vsn = ps.vs_data.vs_num
        phase += vsn
        maxvs = max(maxvs, vsn)
        mass += rho * volume
        energy += ps.w[end] * volume
        minrho = min(minrho, rho)
        maxrho = max(maxrho, rho)
        maxspeed = max(maxspeed, speed)
    end
    return rows, i, phase, maxvs, mass, energy, minrho, maxrho, maxspeed
end

function run_mode(mode::AbstractString)
    spec = mode_spec(mode)
    print_header(spec)

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    nrank = MPI.Comm_size(comm)
    label = envs("VJE_LABEL", spec.mode)
    rundir = joinpath(OUTDIR, label)
    anim_path = joinpath(rundir, "anim")

    config = build_config(spec)
    p4est = nothing
    ka = nothing
    try
        p4est, ka = initialize(config; prerefine_steps = spec.prerefine_steps,
                               prerefine_reinit_ic = true)
        rank == 0 && mkpath(anim_path)
        MPI.Barrier(comm)

        MPI.Barrier(comm)
        t0 = MPI.Wtime()
        solve!(p4est, ka;
            ps_interval = PS_INTERVAL,
            vs_interval = VS_INTERVAL,
            partition_interval = PARTITION_INTERVAL,
            vs_balance = true,
            break_on_convergence = false,
            listen_for_save = false,
            animation = ANIM_DT > 0,
            anim_path = anim_path,
            status_check = true,
            progress = true,
            max_steps = MAX_STEP)
        MPI.Barrier(comm)
        wall = MPI.Wtime() - t0

        rows, loc_ps, loc_phase, loc_maxvs, loc_mass, loc_energy,
            loc_minrho, loc_maxrho, loc_maxspeed = collect_field(ka)

        g_ps = MPI.Reduce(loc_ps, +, comm)
        g_phase = MPI.Reduce(loc_phase, +, comm)
        g_maxvs = MPI.Reduce(loc_maxvs, max, comm)
        g_mass = MPI.Reduce(loc_mass, +, comm)
        g_energy = MPI.Reduce(loc_energy, +, comm)
        g_minrho = MPI.Reduce(loc_minrho, min, comm)
        g_maxrho = MPI.Reduce(loc_maxrho, max, comm)
        g_maxspeed = MPI.Reduce(loc_maxspeed, max, comm)
        data_b = Float64(Base.summarysize(ka.kdata) + Int(KitAMR.p4est_memory_used(p4est)))
        g_data_s = MPI.Reduce(data_b, +, comm)
        g_data_m = MPI.Reduce(data_b, max, comm)
        rss = Float64(Sys.maxrss())
        g_rss_s = MPI.Reduce(rss, +, comm)
        g_rss_m = MPI.Reduce(rss, max, comm)
        n_steps = ka.kinfo.status.step
        sim_time = ka.kinfo.status.sim_time

        rank == 0 && mkpath(rundir)
        MPI.Barrier(comm)
        jldsave(joinpath(rundir, "field_rank$(rank).jld2"); rows = rows)

        if rank == 0
            metrics = (
                mode = spec.mode,
                label = label,
                dimension = DIM,
                np = nrank,
                ps_interval = string(PS_INTERVAL),
                vs_interval = string(VS_INTERVAL),
                partition_interval = string(PARTITION_INTERVAL),
                final_ps_interval_cached = ka.kinfo.status.ps_interval_cached,
                final_vs_interval_cached = ka.kinfo.status.vs_interval_cached,
                final_amr_transport_rate = ka.kinfo.status.amr_transport_rate,
                ps_base = PS_BASE,
                ps_base_x = PS_BASE_X,
                ps_base_y = PS_BASE_Y,
                ps_base_z = DIM == 3 ? PS_BASE_Z : 1,
                ps_lmax = PS_LMAX,
                vs_base = VS_BASE,
                vs_lmax = VS_LMAX,
                min_vs_per_dim = MIN_VS_PER_DIM,
                ps_dynamic = spec.ps_dynamic,
                vs_dynamic = spec.vs_dynamic,
                ps_trees_x = spec.ps_trees_x,
                ps_trees_y = spec.ps_trees_y,
                ps_trees_z = DIM == 3 ? spec.ps_trees_z : 1,
                ps_root_trees = DIM == 2 ? spec.ps_trees_x * spec.ps_trees_y :
                                spec.ps_trees_x * spec.ps_trees_y * spec.ps_trees_z,
                ps_amr_level = spec.ps_amr_level,
                ps_leaves_per_root = 2^(DIM * spec.ps_amr_level),
                both_off_roots_per_rank = BOTH_OFF_ROOTS_PER_RANK,
                both_off_min_roots = BOTH_OFF_MIN_ROOTS,
                vs_trees = spec.vs_trees,
                rho_jet = RHO_JET,
                u_jet = U_JET,
                v_jet = V_JET,
                w_jet = DIM == 3 ? W_JET : 0.0,
                lambda_jet = LAMBDA_JET,
                rho_vac = RHO_VAC,
                lambda_vac = LAMBDA_VAC,
                kn = KN,
                jet_radius = JET_RADIUS,
                jet_edge_width = JET_EDGE_W,
                theta_min = THETA_MIN,
                theta_max = THETA_MAX,
                sigma_min = SIGMA_MIN,
                sigma_max = SIGMA_MAX,
                vrange = VRANGE,
                domain_travel_time = DOMAIN_TRAVEL_TIME,
                domain_l = DOMAIN_L,
                ps_adapt_criterion = "temperature_lohner",
                geometry = GEOM,
                ref_ps_x = REF_PS_X,
                ref_ps_y = REF_PS_Y,
                ref_ps_z = DIM == 3 ? REF_PS_Z : 1,
                ref_vs_per_dim = REF_VS_PER_DIM,
                max_sim_time = MAX_T,
                anim_dt = ANIM_DT,
                n_steps = n_steps,
                sim_time = sim_time,
                wall_s = wall,
                core_h = wall * nrank / 3600,
                total_ps = g_ps,
                total_phase = g_phase,
                max_vs = g_maxvs,
                data_sum_GB = g_data_s / 2^30,
                data_max_GB = g_data_m / 2^30,
                mem_sum_GB = g_rss_s / 2^30,
                mem_max_GB = g_rss_m / 2^30,
                mass = g_mass,
                energy = g_energy,
                minrho = g_minrho,
                maxrho = g_maxrho,
                maxspeed = g_maxspeed,
            )
            jldsave(joinpath(rundir, "metrics.jld2"); metrics = metrics)
            open(joinpath(rundir, "summary.tsv"), "w") do io
                @printf(io, "label\tmode\tnp\tsteps\twall_s\tcore_h\tPS\tphase\tmaxVS\tdata_sum_GB\tdata_max_GB\tmem_sum_GB\tmem_max_GB\tmass\tenergy\n")
                @printf(io, "%s\t%s\t%d\t%d\t%.3f\t%.6f\t%d\t%d\t%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.8e\t%.8e\n",
                        label, spec.mode, nrank, n_steps, wall, wall * nrank / 3600,
                        g_ps, g_phase, g_maxvs, g_data_s / 2^30, g_data_m / 2^30,
                        g_rss_s / 2^30, g_rss_m / 2^30, g_mass, g_energy)
            end
            println("-"^78)
            @printf("[%s] DONE %s : steps=%d wall=%.3fs core-h=%.6f\n",
                    LOG_TAG, label, n_steps, wall, wall * nrank / 3600)
            @printf("[%s]   PS=%d phase=%d maxVS=%d\n", LOG_TAG, g_ps, g_phase, g_maxvs)
            @printf("[%s]   data sum/max=%.3f/%.3f GB; RSS sum/max=%.3f/%.3f GB\n",
                    LOG_TAG,
                    g_data_s / 2^30, g_data_m / 2^30, g_rss_s / 2^30, g_rss_m / 2^30)
            @printf("[%s]   mass=%.8e energy=%.8e rho=[%.3e, %.3e] max|u|=%.3f\n",
                    LOG_TAG,
                    g_mass, g_energy, g_minrho, g_maxrho, g_maxspeed)
            println("-"^78)
        end

        vtk_dir = startswith(OUTDIR, "./") ? OUTDIR[3:end] : OUTDIR
        save_result(p4est, ka; dir_path = joinpath(vtk_dir, label, "result"))
        return nothing
    finally
        if p4est !== nothing && ka !== nothing
            KitAMR.finalize!(p4est, ka)
        end
    end
end

end # module
