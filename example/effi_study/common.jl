module R2DVJEffi

using KitAMR
using MPI
using JLD2
using Printf

const MODES = ("both_on", "ps_only", "vs_only", "both_off")

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

# Fixed vacuum-jet benchmark geometry and states, distilled from
# dev/example/r2d/r2d_vacuum_jet.jl.  The inactive branches in that script
# (bottom jet, fixed-vacuum left boundary, seed packet, and multi-probe VS
# output) are intentionally omitted here.
const GEOM = [0.0, 4.0, -2.0, 2.0]
const XC0 = GEOM[1]
const XC1 = GEOM[2]
const YC0 = GEOM[3]
const YC1 = GEOM[4]
const YMID = 0.5 * (YC0 + YC1)

const PS_BASE = envi("VJE_PS_BASE", 16)
const PS_LMAX = envi("VJE_PS_LMAX", 5)
const VS_BASE = envi("VJE_VS_BASE", 8)
const CFL = envf("VJE_CFL", 0.4)
const MAX_T = envf("VJE_MAX_SIM_TIME", 0.2)
const ANIM_DT = envf("VJE_ANIM_DT", 0.05)
const MAX_STEP = haskey(OPTS, "VJE_MAX_STEPS") ? envi("VJE_MAX_STEPS", 0) : typemax(Int)
const OUTDIR = envs("VJE_OUTDIR", joinpath(@__DIR__, "out"))
const PS_INTERVAL = env_interval("VJE_PS_INTERVAL", "auto")
const VS_INTERVAL = env_interval("VJE_VS_INTERVAL", "auto")
const PARTITION_INTERVAL = env_interval("VJE_PARTITION_INTERVAL", "auto")
const AUTO_PS_TRAVEL_FRACTION = envf("VJE_AUTO_PS_TRAVEL_FRACTION", 2.0)

const RHO_JET = envf("VJE_RHO_JET", 1.0)
const U_JET = envf("VJE_U_JET", 6.0)
const V_JET = envf("VJE_V_JET", 0.0)
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
const U_MAX = hypot(U_JET, V_JET)
const THETA_MAX = THETA_STAT + 0.5 * U_MAX^2 * (GAMMA - 1.0) / GAMMA
const SIGMA_MAX = sqrt(THETA_MAX)
const VRANGE = haskey(OPTS, "VJE_VRANGE") ? envf("VJE_VRANGE", 0.0) : U_MAX + 5.0 * SIGMA_MAX
const VS_LMAX = haskey(OPTS, "VJE_VS_LMAX") ? envi("VJE_VS_LMAX", 0) : 4

const REF_PS_X = PS_BASE * 2^PS_LMAX
const REF_PS_Y = PS_BASE * 2^PS_LMAX
const REF_VS_PER_DIM = VS_BASE * 2^VS_LMAX
const REF_PS_CELLS = REF_PS_X * REF_PS_Y
const REF_VS_CELLS = REF_VS_PER_DIM^2
const DENSE_PHASE_EST = Float64(REF_PS_CELLS) * Float64(REF_VS_CELLS)
const ALLOW_DENSE = envb("VJE_ALLOW_DENSE", false)
const MAX_UNIFORM_PHASE = envf("VJE_MAX_UNIFORM_PHASE", 3.0e8)

function check_inputs()
    PS_BASE > 0 || error("VJE_PS_BASE must be positive")
    PS_LMAX >= 0 || error("VJE_PS_LMAX must be non-negative")
    VS_BASE > 0 || error("VJE_VS_BASE must be positive")
    VS_LMAX >= 0 || error("VJE_VS_LMAX must be non-negative")
    RHO_JET > 0 || error("VJE_RHO_JET must be positive")
    RHO_VAC > 0 || error("VJE_RHO_VAC must be positive")
    LAMBDA_JET > 0 || error("VJE_LAMBDA_JET must be positive")
    LAMBDA_VAC > 0 || error("VJE_LAMBDA_VAC must be positive")
    JET_RADIUS > 0 || error("VJE_JET_RADIUS must be positive")
    JET_EDGE_W > 0 || error("VJE_JET_EDGE_WIDTH must be positive")
    AUTO_PS_TRAVEL_FRACTION > 0 || error("VJE_AUTO_PS_TRAVEL_FRACTION must be positive")
    (YMID - JET_RADIUS > YC0 && YMID + JET_RADIUS < YC1) ||
        error("jet radius must fit inside y-domain")
    return nothing
end

check_inputs()

function smooth_aperture_profile(s, center, radius, edge_width)
    return clamp(0.5 * (1.0 - tanh((abs(s - center) - radius) / edge_width)), 0.0, 1.0)
end

jet_profile_y(y) = smooth_aperture_profile(y, YMID, JET_RADIUS, JET_EDGE_W)

function jet_primitive(a)
    rho = RHO_VAC + a * (RHO_JET - RHO_VAC)
    ux = a * U_JET
    uy = a * V_JET
    theta = THETA_VAC + a * (THETA_JET - THETA_VAC)
    lambda = 1.0 / (2.0 * theta)
    return [rho, ux, uy, lambda]
end

left_jet_primitive(a) = jet_primitive(a)
vacuum_jet_init(midpoint, kinfo) = left_jet_primitive(0.0)
left_aperture_weight(midpoint) = jet_profile_y(midpoint[2])
left_composite_weights(midpoint) = (a = left_aperture_weight(midpoint); (a, 1.0 - a))
const LEFT_JET_BC_PRIM = left_jet_primitive(1.0)

no_vs(ps_data, ka) = false
no_vs(; ps_data, ka) = (0, false)

struct ModeSpec
    mode::String
    ps_dynamic::Bool
    vs_dynamic::Bool
    ps_trees_x::Int
    ps_trees_y::Int
    vs_trees::Int
    ps_amr_level::Int
    vs_amr_level::Int
    prerefine_steps::Int
end

function mode_spec(mode::AbstractString)
    mode in MODES || error("mode must be one of $(join(MODES, ", ")), got $mode")
    ps_dynamic = mode in ("both_on", "ps_only")
    vs_dynamic = mode in ("both_on", "vs_only")
    ps_trees = ps_dynamic ? PS_BASE : PS_BASE * 2^PS_LMAX
    vs_trees = vs_dynamic ? VS_BASE : VS_BASE * 2^VS_LMAX
    ps_amr = ps_dynamic ? PS_LMAX : 0
    vs_amr = vs_dynamic ? VS_LMAX : 0
    prerefine = ps_dynamic ? PS_LMAX : 0
    if mode != "both_on" && DENSE_PHASE_EST > MAX_UNIFORM_PHASE && !ALLOW_DENSE
        error("mode=$mode estimates $(DENSE_PHASE_EST) dense phase cells; lower VJE_PS_LMAX/VJE_VS_LMAX or set VJE_ALLOW_DENSE=1 after checking memory")
    end
    return ModeSpec(String(mode), ps_dynamic, vs_dynamic, ps_trees, ps_trees,
                    vs_trees, ps_amr, vs_amr, prerefine)
end

function build_config(spec::ModeSpec)
    solver = Solver(;
        DIM = 2,
        NDF = 2,
        CFL = CFL,
        AMR_PS_MAXLEVEL = spec.ps_amr_level,
        AMR_VS_MAXLEVEL = spec.vs_amr_level,
        AMR_PS_DYNAMIC = spec.ps_dynamic,
        AMR_PS_THRES = 0.3,
        AMR_PS_SMOOTH = 0.5,
        AMR_VS_DYNAMIC = spec.vs_dynamic,
        AUTO_AMR_PS_TRAVEL_FRACTION = AUTO_PS_TRAVEL_FRACTION,
        flux = CAIDVM,
        time_marching = CIP_Marching,
        max_sim_time = MAX_T,
    )

    gas = Gas(; K = 1.0, Kn = KN, ω = 0.81, ωᵣ = 0.81)
    output = Output(solver;
        vtk_celltype = [Pixel, Triangle],
        vs_vtk_celltype = Pixel,
        anim_dt = ANIM_DT,
        vs_output_criterion = no_vs,
    )
    left_domain = Domain(Composite, 1, CompositeBC(left_composite_weights,
        Domain(SuperSonicInflow, 1, LEFT_JET_BC_PRIM),
        Domain(UniformOutflow, 1)))

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
        user_defined = UDF(),
    )
end

function print_header(spec::ModeSpec)
    rank = MPI.Comm_rank(MPI.COMM_WORLD)
    rank == 0 || return nothing
    println("="^78)
    @printf("[r2d_vj_effi] mode=%s np=%d\n", spec.mode, MPI.Comm_size(MPI.COMM_WORLD))
    @printf("[r2d_vj_effi] geometry=[%.3f, %.3f] x [%.3f, %.3f]\n",
            GEOM[1], GEOM[2], GEOM[3], GEOM[4])
    @printf("[r2d_vj_effi] PS trees=(%d,%d) AMR_L=%d dynamic=%s finest=%d x %d\n",
            spec.ps_trees_x, spec.ps_trees_y, spec.ps_amr_level,
            string(spec.ps_dynamic), REF_PS_X, REF_PS_Y)
    @printf("[r2d_vj_effi] VS trees=%d AMR_L=%d dynamic=%s finest=%d/dim\n",
            spec.vs_trees, spec.vs_amr_level, string(spec.vs_dynamic), REF_VS_PER_DIM)
    @printf("[r2d_vj_effi] max_t=%.4g anim_dt=%.4g vrange=+/-%.4g Kn=%.4g\n",
            MAX_T, ANIM_DT, VRANGE, KN)
    @printf("[r2d_vj_effi] AMR intervals: PS=%s VS=%s partition=%s vs_balance=true\n",
            string(PS_INTERVAL), string(VS_INTERVAL), string(PARTITION_INTERVAL))
    @printf("[r2d_vj_effi] outdir=%s\n", OUTDIR)
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

    rows = Matrix{Float64}(undef, nloc, 8) # x y dx dy rho ux uy theta
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
        area = prod(ps.ds)
        rho = ps.prim[1]
        theta = 1.0 / (2.0 * ps.prim[end])
        speed = hypot(ps.prim[2], ps.prim[3])
        rows[i, 1] = ps.midpoint[1]
        rows[i, 2] = ps.midpoint[2]
        rows[i, 3] = ps.ds[1]
        rows[i, 4] = ps.ds[2]
        rows[i, 5] = rho
        rows[i, 6] = ps.prim[2]
        rows[i, 7] = ps.prim[3]
        rows[i, 8] = theta
        vsn = ps.vs_data.vs_num
        phase += vsn
        maxvs = max(maxvs, vsn)
        mass += rho * area
        energy += ps.w[end] * area
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
                np = nrank,
                ps_interval = string(PS_INTERVAL),
                vs_interval = string(VS_INTERVAL),
                partition_interval = string(PARTITION_INTERVAL),
                final_ps_interval_cached = ka.kinfo.status.ps_interval_cached,
                final_vs_interval_cached = ka.kinfo.status.vs_interval_cached,
                final_amr_transport_rate = ka.kinfo.status.amr_transport_rate,
                ps_base = PS_BASE,
                ps_lmax = PS_LMAX,
                vs_base = VS_BASE,
                vs_lmax = VS_LMAX,
                ps_dynamic = spec.ps_dynamic,
                vs_dynamic = spec.vs_dynamic,
                ps_trees_x = spec.ps_trees_x,
                ps_trees_y = spec.ps_trees_y,
                vs_trees = spec.vs_trees,
                rho_jet = RHO_JET,
                u_jet = U_JET,
                v_jet = V_JET,
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
                geometry = GEOM,
                ref_ps_x = REF_PS_X,
                ref_ps_y = REF_PS_Y,
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
            @printf("[r2d_vj_effi] DONE %s : steps=%d wall=%.3fs core-h=%.6f\n",
                    label, n_steps, wall, wall * nrank / 3600)
            @printf("[r2d_vj_effi]   PS=%d phase=%d maxVS=%d\n", g_ps, g_phase, g_maxvs)
            @printf("[r2d_vj_effi]   data sum/max=%.3f/%.3f GB; RSS sum/max=%.3f/%.3f GB\n",
                    g_data_s / 2^30, g_data_m / 2^30, g_rss_s / 2^30, g_rss_m / 2^30)
            @printf("[r2d_vj_effi]   mass=%.8e energy=%.8e rho=[%.3e, %.3e] max|u|=%.3f\n",
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
