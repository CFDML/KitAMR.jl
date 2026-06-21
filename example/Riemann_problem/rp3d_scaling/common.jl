module RP3DScaling

using KitAMR
using MPI
using Printf

const NODE_MEMORY_GB = 256.0
const JULIA_BASELINE_GB_PER_RANK = 0.8
const BYTES_PER_PHASE_LOW = 250.0
const BYTES_PER_PHASE_HIGH = 450.0

struct CaseSpec
    label::String
    scaling::Symbol
    expected_cores::Int
    trees_num::NTuple{3,Int}
    ps_maxlevel::Int
    ps_dynamic_maxlevel::Int
    use_uniform_static_refine::Bool
    uniform_refine_level::Int
    amr_ps_dynamic::Bool
    prerefine_steps::Int
    force_partition::Bool
    vs_trees_num::NTuple{3,Int}
    vs_maxlevel::Int
    amr_vs_dynamic::Bool
    warmup_steps::Int
    timed_steps::Int
    cfl::Float64
    amr_ps_threshold::Float64
    quadrature_bound::Float64
end

env_int(name::AbstractString, default::Integer) =
    something(tryparse(Int, get(ENV, name, "")), Int(default))

env_float(name::AbstractString, default::Real) =
    something(tryparse(Float64, get(ENV, name, "")), Float64(default))

env_bool(name::AbstractString, default::Bool) =
    lowercase(get(ENV, name, default ? "1" : "0")) in ("1", "true", "yes", "on")

env_optional_float(name::AbstractString) =
    haskey(ENV, name) ? tryparse(Float64, ENV[name]) : nothing

const DEFAULT_BULK_SPEED_BOUND = 0.7276
const DEFAULT_SIGMA_MULTIPLIER = 5.0
const DEFAULT_TMAX_ESTIMATE = 2.0
const DEFAULT_SETUP_TIME = 0.01
const STRONG_VS_NUM_ESTIMATE = (
    mean = 6185.5,
)

function estimated_quadrature_bound()
    bulk = env_float("RP3D_BULK_SPEED_BOUND", DEFAULT_BULK_SPEED_BOUND)
    nsigma = env_float("RP3D_SIGMA_MULTIPLIER", DEFAULT_SIGMA_MULTIPLIER)
    tmax = env_float("RP3D_TMAX_ESTIMATE", DEFAULT_TMAX_ESTIMATE)
    return bulk + nsigma * sqrt(tmax)
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

lambda_from_pressure(rho, p) = 0.5 * rho / p

state_from_rho_u_v_w_p(rho, u, v, w, p) =
    [rho, u, v, w, lambda_from_pressure(rho, p)]

function riemann_3d_init(midpoint, kinfo)
    x = midpoint[1]
    y = midpoint[2]
    z = midpoint[3]

    if z >= 0.0
        if x >= 0.0 && y >= 0.0
            return state_from_rho_u_v_w_p(0.5313, 0.0, 0.0, 0.0, 0.4)
        elseif x < 0.0 && y >= 0.0
            return state_from_rho_u_v_w_p(1.0, 0.7276, 0.0, 0.0, 1.0)
        elseif x < 0.0 && y < 0.0
            return state_from_rho_u_v_w_p(0.8, 0.0, 0.0, 0.0, 1.0)
        else
            return state_from_rho_u_v_w_p(1.0, 0.0, 0.7276, 0.0, 1.0)
        end
    else
        if x >= 0.0 && y >= 0.0
            return state_from_rho_u_v_w_p(1.0, 0.0, 0.0, 0.7276, 1.0)
        elseif x < 0.0 && y >= 0.0
            return state_from_rho_u_v_w_p(0.8, 0.7276, 0.0, 0.7276, 1.0)
        elseif x < 0.0 && y < 0.0
            return state_from_rho_u_v_w_p(0.5313, 0.0, 0.0, 0.7276, 0.4)
        else
            return state_from_rho_u_v_w_p(0.8, 0.0, 0.7276, 0.7276, 1.0)
        end
    end
end

function strong_spec(cores::Integer)
    return CaseSpec(
        "rp3d_strong_$(cores)c",
        :strong,
        Int(cores),
        (16, 16, 16),
        4,
        4,
        false,
        0,
        true,
        4,
        true,
        (4, 4, 4),
        4,
        false,
        8,
        20,
        0.4,
        0.3,
        estimated_quadrature_bound(),
    )
end

const WEAK_ROOT_TREES = Dict(
    64 => 6,
    128 => 7,
    256 => 9,
    512 => 11,
    1024 => 14,
    2048 => 18,
)

function weak_spec(cores::Integer)
    c = Int(cores)
    haskey(WEAK_ROOT_TREES, c) ||
        error("No weak-scaling grid preset for $c cores.")
    root = WEAK_ROOT_TREES[c]
    return CaseSpec(
        "rp3d_weak_$(c)c",
        :weak,
        c,
        (root, root, root),
        3,
        3,
        true,
        3,
        false,
        0,
        true,
        (16, 16, 16),
        0,
        false,
        5,
        40,
        0.4,
        0.3,
        estimated_quadrature_bound(),
    )
end

function quick_probe_spec()
    return CaseSpec(
        "rp3d_quick_probe_2048c",
        :probe,
        2048,
        (13, 13, 13),
        3,
        3,
        false,
        0,
        true,
        3,
        true,
        (4, 4, 4),
        4,
        false,
        1,
        3,
        0.4,
        0.3,
        estimated_quadrature_bound(),
    )
end

function precompile_strong_spec()
    return CaseSpec(
        "rp3d_precompile_strong",
        :precompile,
        1,
        (2, 2, 2),
        2,
        2,
        false,
        0,
        true,
        1,
        true,
        (8, 8, 8),
        0,
        false,
        1,
        1,
        0.4,
        0.3,
        5.0,
    )
end

function precompile_weak_spec()
    return CaseSpec(
        "rp3d_precompile_weak",
        :precompile,
        1,
        (2, 2, 2),
        1,
        1,
        true,
        1,
        false,
        0,
        true,
        (8, 8, 8),
        0,
        false,
        1,
        1,
        0.4,
        0.3,
        5.0,
    )
end

function with_env_overrides(spec::CaseSpec)
    warmup_steps = env_int("RP3D_WARMUP_STEPS", spec.warmup_steps)
    timed_steps = env_int("RP3D_TIMED_STEPS", spec.timed_steps)
    qbound = env_float("RP3D_QUADRATURE_BOUND", spec.quadrature_bound)
    return CaseSpec(
        spec.label,
        spec.scaling,
        spec.expected_cores,
        spec.trees_num,
        spec.ps_maxlevel,
        spec.ps_dynamic_maxlevel,
        spec.use_uniform_static_refine,
        spec.uniform_refine_level,
        spec.amr_ps_dynamic,
        spec.prerefine_steps,
        spec.force_partition,
        spec.vs_trees_num,
        spec.vs_maxlevel,
        spec.amr_vs_dynamic,
        warmup_steps,
        timed_steps,
        spec.cfl,
        spec.amr_ps_threshold,
        qbound,
    )
end

function uniform_static_refine_flag(target_level::Integer)
    return (midpoint, ds, kinfo, level) -> level < target_level
end

function build_config(spec::CaseSpec)
    solver = Solver(;
        DIM = 3,
        NDF = 1,
        CFL = spec.cfl,
        AMR_PS_MAXLEVEL = spec.ps_maxlevel,
        AMR_PS_DYNAMIC_MAXLEVEL = spec.ps_dynamic_maxlevel,
        AMR_VS_MAXLEVEL = spec.vs_maxlevel,
        AMR_PS_DYNAMIC = spec.amr_ps_dynamic,
        AMR_PS_THRES = spec.amr_ps_threshold,
        AMR_VS_DYNAMIC = spec.amr_vs_dynamic,
        flux = CAIDVM,
        time_marching = CIP_Marching,
        ST_CHECK_INTERVAL = typemax(Int),
        max_sim_time = Inf,
    )

    gas = Gas(;
        K = 0.0,
        Kn = 0.001,
        ω = 0.81,
        ωᵣ = 0.81,
    )

    output = Output(solver; anim_dt = 0.0)
    udf = spec.use_uniform_static_refine ?
        UDF(; static_ps_refine_flag = uniform_static_refine_flag(spec.uniform_refine_level)) :
        UDF()
    q = spec.quadrature_bound

    return Configure(solver;
        geometry = [-0.5, 0.5, -0.5, 0.5, -0.5, 0.5],
        trees_num = collect(spec.trees_num),
        quadrature = [-q, q, -q, q, -q, q],
        vs_trees_num = collect(spec.vs_trees_num),
        IC = PCoordFn(riemann_3d_init),
        domain = [
            Domain(UniformOutflow, 1), Domain(UniformOutflow, 2),
            Domain(UniformOutflow, 3), Domain(UniformOutflow, 4),
            Domain(UniformOutflow, 5), Domain(UniformOutflow, 6),
        ],
        output = output,
        gas = gas,
        user_defined = udf,
    )
end

function local_stats(ka, p4est)
    comm = MPI.COMM_WORLD
    local_ps = 0
    local_phase = 0
    local_max_vs = 0
    local_max_temperature = 0.0
    local_max_component_velocity = 0.0
    for tree in ka.kdata.field.trees.data
        for ps_data in tree
            ps_data isa InsideSolidData && continue
            vs_num = ps_data.vs_data.vs_num
            local_ps += 1
            local_phase += vs_num
            local_max_vs = max(local_max_vs, vs_num)
            prim = ps_data.prim
            local_max_temperature = max(local_max_temperature, 1.0 / (2.0 * prim[end]))
            local_max_component_velocity = max(local_max_component_velocity,
                maximum(abs, @view(prim[2:4])))
        end
    end

    total_ps = MPI.Allreduce(local_ps, +, comm)
    total_phase = MPI.Allreduce(local_phase, +, comm)
    max_vs = MPI.Allreduce(local_max_vs, MPI.MAX, comm)
    max_temperature = MPI.Allreduce(local_max_temperature, MPI.MAX, comm)
    max_component_velocity = MPI.Allreduce(local_max_component_velocity, MPI.MAX, comm)
    min_local_phase = MPI.Allreduce(local_phase, MPI.MIN, comm)
    max_local_phase = MPI.Allreduce(local_phase, MPI.MAX, comm)
    sum_weight = MPI.Allreduce(Float64(local_phase), +, comm)
    max_weight = MPI.Allreduce(Float64(local_phase), MPI.MAX, comm)
    nranks = MPI.Comm_size(comm)
    imbalance = sum_weight > 0.0 ? max_weight / (sum_weight / nranks) - 1.0 : 0.0
    return (
        total_ps = total_ps,
        global_quads = total_ps,
        total_phase = total_phase,
        max_vs = max_vs,
        max_temperature = max_temperature,
        max_component_velocity = max_component_velocity,
        recommended_q_5sigma = max_component_velocity + DEFAULT_SIGMA_MULTIPLIER * sqrt(max_temperature),
        min_local_phase = min_local_phase,
        max_local_phase = max_local_phase,
        imbalance = imbalance,
    )
end

function estimate_memory(total_phase::Integer, nranks::Integer)
    data_low_gb = total_phase * BYTES_PER_PHASE_LOW / 1e9
    data_high_gb = total_phase * BYTES_PER_PHASE_HIGH / 1e9
    baseline_gb = nranks * JULIA_BASELINE_GB_PER_RANK
    total_low_gb = data_low_gb + baseline_gb
    total_high_gb = data_high_gb + baseline_gb
    nodes_low = ceil(Int, total_low_gb / NODE_MEMORY_GB)
    nodes_high = ceil(Int, total_high_gb / NODE_MEMORY_GB)
    return (
        data_low_gb = data_low_gb,
        data_high_gb = data_high_gb,
        total_low_gb = total_low_gb,
        total_high_gb = total_high_gb,
        nodes_low = nodes_low,
        nodes_high = nodes_high,
    )
end

function print_case_header(spec::CaseSpec)
    rank = MPI.Comm_rank(MPI.COMM_WORLD)
    rank == 0 || return nothing
    final_n = spec.trees_num[1] * 2^spec.ps_maxlevel
    root_trees = prod(spec.trees_num)
    ps_mode = spec.use_uniform_static_refine ?
        @sprintf("uniform_static(level=%d)", spec.uniform_refine_level) :
        "initial_amr"
    vs_mode = spec.vs_maxlevel == 0 && !spec.amr_vs_dynamic ?
        "uniform" :
        (spec.amr_vs_dynamic ? "dynamic_amr" : "initial_maxwellian_amr")
    @printf("\n=== %s ===\n", spec.label)
    @printf("scaling=%s expected_cores=%d mpi_ranks=%d\n",
        String(spec.scaling), spec.expected_cores, MPI.Comm_size(MPI.COMM_WORLD))
    @printf("physical_grid=%s AMR_PS_DYNAMIC=%s velocity_grid=%s\n",
        ps_mode, string(spec.amr_ps_dynamic), vs_mode)
    @printf("root_trees=%d\n", root_trees)
    if root_trees < MPI.Comm_size(MPI.COMM_WORLD)
        @warn "root_trees is smaller than mpi_ranks; AMR/partition can become rank-count dependent" root_trees mpi_ranks = MPI.Comm_size(MPI.COMM_WORLD)
    end
    @printf("trees_num=[%d,%d,%d] AMR_PS_MAXLEVEL=%d finest_equiv_N=%d\n",
        spec.trees_num..., spec.ps_maxlevel, final_n)
    @printf("vs_trees_num=[%d,%d,%d] AMR_VS_MAXLEVEL=%d AMR_VS_DYNAMIC=%s warmup=%d timed=%d\n",
        spec.vs_trees_num..., spec.vs_maxlevel, string(spec.amr_vs_dynamic),
        spec.warmup_steps, spec.timed_steps)
    @printf("quadrature=[-%.4g, %.4g]^3 (set RP3D_TMAX_ESTIMATE or RP3D_QUADRATURE_BOUND to override)\n",
        spec.quadrature_bound, spec.quadrature_bound)
    return nothing
end

function print_stats(label::AbstractString, stats, nranks::Integer)
    rank = MPI.Comm_rank(MPI.COMM_WORLD)
    rank == 0 || return nothing
    mem = estimate_memory(stats.total_phase, nranks)
    @printf("[%s] physical=%d phase=%d max_vs=%d imbalance=%.3f local_phase=[%d,%d]\n",
        label, stats.global_quads, stats.total_phase, stats.max_vs, stats.imbalance,
        stats.min_local_phase, stats.max_local_phase)
    @printf("[%s] max_T=%.6g max_component_U=%.6g recommended_q_5sigma=%.6g\n",
        label, stats.max_temperature, stats.max_component_velocity, stats.recommended_q_5sigma)
    @printf("[%s] estimated data=%.2f-%.2f TB, incl Julia baseline=%.2f-%.2f TB, min nodes=%d-%d\n",
        label,
        mem.data_low_gb / 1000.0,
        mem.data_high_gb / 1000.0,
        mem.total_low_gb / 1000.0,
        mem.total_high_gb / 1000.0,
        mem.nodes_low,
        mem.nodes_high)
    return nothing
end

function timed_barrier_elapsed(f::Function)
    comm = MPI.COMM_WORLD
    MPI.Barrier(comm)
    t0 = time()
    f()
    MPI.Barrier(comm)
    elapsed = time() - t0
    return MPI.Allreduce(elapsed, MPI.MAX, comm)
end

function manual_gc_before_timing!(label::AbstractString)
    env_bool("RP3D_GC_BEFORE_TIMING", true) || return 0.0
    passes = env_int("RP3D_GC_PASSES", 1)
    passes > 0 || return 0.0
    elapsed = timed_barrier_elapsed() do
        for _ in 1:passes
            GC.gc(true)
        end
    end
    rank = MPI.Comm_rank(MPI.COMM_WORLD)
    rank == 0 &&
        @printf("[manual_gc] label=%s passes=%d elapsed=%.3f s\n", label, passes, elapsed)
    return elapsed
end

function march_steps!(p4est, ka, steps::Integer)
    for _ in 1:steps
        limit_dt!(ka)
        slope!(p4est, ka)
        flux!(p4est, ka)
        iterate!(p4est, ka)
    end
    return nothing
end

limit_dt!(ka) = KitAMR.limit_Δt!(ka)

function force_partition!(p4est, ka)
    ps_partition!(p4est, ka)
    amr_recover!(p4est, ka; topology_changed = true, velocity_changed = true)
    return nothing
end

function run_benchmark(input_spec::CaseSpec; allow_wrong_cores::Bool = false)
    spec = with_env_overrides(input_spec)
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    nranks = MPI.Comm_size(comm)
    if rank == 0 && nranks != spec.expected_cores
        msg = "script preset expects $(spec.expected_cores) MPI ranks, got $nranks"
        allow_wrong_cores ? @warn(msg) : @warn(msg * " (continuing anyway)")
    end

    print_case_header(spec)
    config = build_config(spec)
    p4est = nothing
    ka = nothing
    try
        init_elapsed = timed_barrier_elapsed() do
            p4est, ka = initialize(config;
                prerefine_steps = spec.prerefine_steps,
                prerefine_reinit_ic = true)
        end
        stats_init = local_stats(ka, p4est)
        print_stats("after_init", stats_init, nranks)

        partition_elapsed = 0.0
        if spec.force_partition
            partition_elapsed = timed_barrier_elapsed() do
                force_partition!(p4est, ka)
            end
            stats_partitioned = local_stats(ka, p4est)
            print_stats("after_partition", stats_partitioned, nranks)
        end

        warmup_elapsed = spec.warmup_steps > 0 ?
            timed_barrier_elapsed(() -> march_steps!(p4est, ka, spec.warmup_steps)) : 0.0

        gc_elapsed = manual_gc_before_timing!("$(spec.label)_before_timed")
        stats_timed = local_stats(ka, p4est)
        elapsed = timed_barrier_elapsed() do
            march_steps!(p4est, ka, spec.timed_steps)
        end
        seconds_per_step = elapsed / max(spec.timed_steps, 1)
        phase_rate = stats_timed.total_phase * spec.timed_steps / max(elapsed, eps())

        if rank == 0
            @printf("[timing] init=%.3f s partition=%.3f s warmup=%.3f s manual_gc=%.3f s timed=%.3f s steps=%d sec_per_step=%.6f phase_updates_per_sec=%.6e\n",
                init_elapsed, partition_elapsed, warmup_elapsed, gc_elapsed, elapsed, spec.timed_steps,
                seconds_per_step, phase_rate)
        end

        return (
            spec = spec,
            stats = stats_timed,
            elapsed = elapsed,
            seconds_per_step = seconds_per_step,
            phase_rate = phase_rate,
            gc_elapsed = gc_elapsed,
        )
    finally
        if p4est !== nothing && ka !== nothing
            finalize!(p4est, ka)
        end
    end
end

function run_initial_amr_benchmark(input_spec::CaseSpec; allow_wrong_cores::Bool = false)
    spec = with_env_overrides(input_spec)
    setup_time = env_float("RP3D_SETUP_TIME", DEFAULT_SETUP_TIME)
    setup_steps_override = env_int("RP3D_SETUP_STEPS", 0)
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    nranks = MPI.Comm_size(comm)
    if rank == 0 && nranks != spec.expected_cores
        msg = "script preset expects $(spec.expected_cores) MPI ranks, got $nranks"
        allow_wrong_cores ? @warn(msg) : @warn(msg * " (continuing anyway)")
    end

    print_case_header(spec)
    config = build_config(spec)
    p4est = nothing
    ka = nothing
    try
        init_elapsed = timed_barrier_elapsed() do
            p4est, ka = initialize(config;
                prerefine_steps = spec.prerefine_steps,
                prerefine_reinit_ic = true)
        end
        print_stats("after_initial_prerefine", local_stats(ka, p4est), nranks)

        setup_steps_done = Ref(0)
        setup_elapsed = timed_barrier_elapsed() do
            if setup_steps_override > 0
                adaptive_march_steps!(p4est, ka, setup_steps_override)
                setup_steps_done[] = setup_steps_override
            elseif setup_time > 0.0
                setup_steps_done[] = adaptive_march_until!(p4est, ka, setup_time)
            end
        end
        print_stats("after_adaptive_setup", local_stats(ka, p4est), nranks)

        partition_elapsed = 0.0
        if spec.force_partition
            partition_elapsed = timed_barrier_elapsed() do
                force_partition!(p4est, ka)
            end
            print_stats("after_final_partition", local_stats(ka, p4est), nranks)
        end

        warmup_elapsed = spec.warmup_steps > 0 ?
            timed_barrier_elapsed(() -> march_steps!(p4est, ka, spec.warmup_steps)) : 0.0

        gc_elapsed = manual_gc_before_timing!("$(spec.label)_before_timed")
        stats_timed = local_stats(ka, p4est)
        elapsed = timed_barrier_elapsed() do
            march_steps!(p4est, ka, spec.timed_steps)
        end
        seconds_per_step = elapsed / max(spec.timed_steps, 1)
        phase_rate = stats_timed.total_phase * spec.timed_steps / max(elapsed, eps())

        if rank == 0
            @printf("[initial_amr_timing] init=%.3f s adaptive_setup=%.3f s setup_time_target=%.6g actual_sim_time=%.6g setup_steps=%d partition=%.3f s warmup=%.3f s manual_gc=%.3f s timed=%.3f s steps=%d sec_per_step=%.6f phase_updates_per_sec=%.6e\n",
                init_elapsed, setup_elapsed, setup_time, ka.kinfo.status.sim_time,
                setup_steps_done[], partition_elapsed, warmup_elapsed, gc_elapsed, elapsed,
                spec.timed_steps, seconds_per_step, phase_rate)
        end

        return (
            spec = spec,
            stats = stats_timed,
            elapsed = elapsed,
            seconds_per_step = seconds_per_step,
            phase_rate = phase_rate,
            gc_elapsed = gc_elapsed,
        )
    finally
        if p4est !== nothing && ka !== nothing
            finalize!(p4est, ka)
        end
    end
end

function adaptive_march_steps!(p4est, ka, steps::Integer;
        ps_interval = 1, vs_interval = 2, partition_interval = 2)
    for _ in 1:steps
        adaptive_mesh_refinement!(p4est, ka;
            ps_interval = ps_interval,
            vs_interval = vs_interval,
            partition_interval = partition_interval)
        limit_dt!(ka)
        slope!(p4est, ka)
        flux!(p4est, ka)
        iterate!(p4est, ka)
    end
    return nothing
end

function adaptive_march_until!(p4est, ka, setup_time::Real;
        max_steps::Integer = typemax(Int),
        ps_interval = 1, vs_interval = 2, partition_interval = 2)
    start_time = ka.kinfo.status.sim_time
    target_time = start_time + setup_time
    steps = 0
    while ka.kinfo.status.sim_time < target_time && steps < max_steps
        adaptive_mesh_refinement!(p4est, ka;
            ps_interval = ps_interval,
            vs_interval = vs_interval,
            partition_interval = partition_interval)
        limit_dt!(ka)
        slope!(p4est, ka)
        flux!(p4est, ka)
        iterate!(p4est, ka)
        steps += 1
    end
    return steps
end

function run_weak_uniform(cores::Integer; allow_wrong_cores::Bool = false)
    return run_benchmark(weak_spec(cores); allow_wrong_cores = allow_wrong_cores)
end

function run_strong_initial_amr(cores::Integer; allow_wrong_cores::Bool = false)
    return run_initial_amr_benchmark(strong_spec(cores); allow_wrong_cores = allow_wrong_cores)
end

function estimate_target_phase_counts()
    strong = strong_spec(2048)
    weak = weak_spec(2048)
    n_eff_strong = strong.trees_num[1] * 2^strong.ps_maxlevel
    coeff_low = env_float("RP3D_STRONG_LEAF_COEFF_LOW", 12.0)
    coeff_high = env_float("RP3D_STRONG_LEAF_COEFF_HIGH", 20.0)
    strong_leaf_low = round(Int, coeff_low * n_eff_strong^2)
    strong_leaf_high = round(Int, coeff_high * n_eff_strong^2)
    strong_vs = strong.vs_maxlevel == 0 ? prod(strong.vs_trees_num) : STRONG_VS_NUM_ESTIMATE.mean
    weak_n = weak.trees_num[1] * 2^weak.ps_maxlevel
    weak_leaf = weak_n^3
    weak_vs = prod(weak.vs_trees_num)
    return (
        strong_leaf_low = strong_leaf_low,
        strong_leaf_high = strong_leaf_high,
        strong_vs_mean = strong_vs,
        strong_phase_low = round(Int, strong_leaf_low * strong_vs),
        strong_phase_high = round(Int, strong_leaf_high * strong_vs),
        weak_leaf = weak_leaf,
        weak_phase = weak_leaf * weak_vs,
    )
end

function print_quick_estimates(result)
    rank = MPI.Comm_rank(MPI.COMM_WORLD)
    rank == 0 || return nothing
    est = estimate_target_phase_counts()
    rate = result.phase_rate
    strong_low_time = est.strong_phase_low / max(rate, eps())
    strong_high_time = est.strong_phase_high / max(rate, eps())
    weak_time = est.weak_phase / max(rate, eps())
    strong_mem_low = estimate_memory(est.strong_phase_low, 2048)
    strong_mem_high = estimate_memory(est.strong_phase_high, 2048)
    weak_mem = estimate_memory(est.weak_phase, 2048)

    @printf("\n=== 2048-core production estimates from quick probe ===\n")
    @printf("strong target physical=%d-%d mean_vs≈%.1f phase=%d-%d estimated_sec_per_step=%.3f-%.3f\n",
        est.strong_leaf_low, est.strong_leaf_high,
        est.strong_vs_mean,
        est.strong_phase_low, est.strong_phase_high,
        strong_low_time, strong_high_time)
    @printf("strong target data=%.2f-%.2f TB, incl baseline=%.2f-%.2f TB, min nodes=%d-%d\n",
        strong_mem_low.data_low_gb / 1000.0,
        strong_mem_high.data_high_gb / 1000.0,
        strong_mem_low.total_low_gb / 1000.0,
        strong_mem_high.total_high_gb / 1000.0,
        strong_mem_low.nodes_low,
        strong_mem_high.nodes_high)
    @printf("weak target physical=%d phase=%d estimated_sec_per_step=%.3f\n",
        est.weak_leaf, est.weak_phase, weak_time)
    @printf("weak target data=%.2f-%.2f TB, incl baseline=%.2f-%.2f TB, min nodes=%d-%d\n",
        weak_mem.data_low_gb / 1000.0,
        weak_mem.data_high_gb / 1000.0,
        weak_mem.total_low_gb / 1000.0,
        weak_mem.total_high_gb / 1000.0,
        weak_mem.nodes_low,
        weak_mem.nodes_high)
    return nothing
end

function run_quick_probe()
    result = run_initial_amr_benchmark(quick_probe_spec(); allow_wrong_cores = true)
    print_quick_estimates(result)
    return result
end

function run_precompile_workload()
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)

    # Small dynamic-AMR path.
    spec = precompile_strong_spec()
    p4est = nothing
    ka = nothing
    try
        config = build_config(spec)
        p4est, ka = initialize(config;
            prerefine_steps = spec.prerefine_steps,
            prerefine_reinit_ic = true)
        adaptive_march_steps!(p4est, ka, 1)
        force_partition!(p4est, ka)
        march_steps!(p4est, ka, 1)
    finally
        if p4est !== nothing && ka !== nothing
            finalize!(p4est, ka)
        end
    end

    # Uniform weak-scaling path.
    p4est = nothing
    ka = nothing
    try
        spec = precompile_weak_spec()
        config = build_config(spec)
        p4est, ka = initialize(config;
            prerefine_steps = spec.prerefine_steps,
            prerefine_reinit_ic = true)
        force_partition!(p4est, ka)
        march_steps!(p4est, ka, 1)
    finally
        if p4est !== nothing && ka !== nothing
            finalize!(p4est, ka)
        end
    end

    rank == 0 && println("RP3D precompile workload complete.")
    return nothing
end

end # module
