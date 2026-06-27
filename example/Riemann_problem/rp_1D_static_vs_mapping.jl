import Dates
using KitAMR
using MPI
using Printf

MPI.Init()

const LEFT_VS_LEVEL = 4
const RIGHT_VS_LEVEL = 3
const SPLIT_X = 0.0
const PHYSICAL_TREES_NUM = [
    parse(Int, get(ENV, "KITAMR_RP1D_NX", "128")),
    parse(Int, get(ENV, "KITAMR_RP1D_NY", "4")),
]
const MAX_SIM_TIME = parse(Float64, get(ENV, "KITAMR_RP1D_MAX_SIM_TIME", "0.2"))
const RUN_ID = get(ENV, "KITAMR_RP1D_RUN_ID",
                   Dates.format(Dates.now(), "yyyymmdd_HHMMSS"))
const OUTPUT_ROOT = joinpath("example", "Riemann_problem", "results",
                             "static_vs_mapping_" * RUN_ID)

function Sod_init(midpoint, kinfo)
    if midpoint[1] < 0.0
        return [1.0, 0.0, 0.0, 1.0]
    else
        return [0.125, 0.0, 0.0, 0.625]
    end
end

function make_config()
    solver = Solver(;
        DIM = 2, NDF = 2,
        CFL = 0.4,
        AMR_PS_MAXLEVEL = 4,
        AMR_VS_MAXLEVEL = LEFT_VS_LEVEL,
        AMR_PS_DYNAMIC = false,
        AMR_PS_THRES = 0.25,
        AMR_VS_DYNAMIC = false,
        flux = CAIDVM,
        time_marching = CIP_Marching,
        ST_CHECK_INTERVAL = 20,
        max_sim_time = MAX_SIM_TIME,
    )
    gas = Gas(;
        K = 1.0,
        Kn = 1e3,
        ω = 0.81,
        ωᵣ = 0.81,
    )
    output = Output(solver)
    udf = UDF()
    return Configure(solver;
        geometry = [-2.0, 2.0, -0.5, 0.5],
        trees_num = PHYSICAL_TREES_NUM,
        quadrature = [-7.0, 7.0, -7.0, 7.0],
        vs_trees_num = [8, 8],
        IC = PCoordFn(Sod_init),
        domain = [
            Domain(UniformOutflow, 1), Domain(UniformOutflow, 2),
            Domain(Period, 3), Domain(Period, 4),
        ],
        output = output,
        gas = gas,
        user_defined = udf,
    )
end

function velocity_cell_size(config)
    return [
        (config.quadrature[2 * d] - config.quadrature[2 * d - 1]) /
        config.vs_trees_num[d]
        for d in 1:2
    ]
end

function cap_velocity_maxlevel!(ka; split_x = SPLIT_X,
                                left_level = LEFT_VS_LEVEL,
                                right_level = RIGHT_VS_LEVEL)
    ds = velocity_cell_size(ka.kinfo.config)
    global_maxlevel = ka.kinfo.config.solver.AMR_VS_MAXLEVEL
    changed_cells = 0

    for tree in ka.kdata.field.trees.data
        for ps_data in tree
            hasproperty(ps_data, :vs_data) || continue
            target_level = ps_data.midpoint[1] < split_x ? left_level : right_level
            vs_data = ps_data.vs_data
            changed = false
            while maximum(vs_data.level) > target_level
                coarsen_ok = vs_data.level .> target_level
                did_coarsen =
                    KitAMR.coarsen_grid_stream!(vs_data, coarsen_ok, ds, global_maxlevel)
                did_coarsen ||
                    error("Could not coarsen velocity grid at physical midpoint ",
                          ps_data.midpoint, " to level ", target_level)
                changed = true
            end
            if changed
                KitAMR.conserved_I_projection!(vs_data, ps_data.w)
                changed_cells += 1
            end
        end
    end

    return MPI.Allreduce(changed_cells, MPI.SUM, MPI.COMM_WORLD)
end

function local_vs_stats(ka)
    left_ps = 0
    right_ps = 0
    left_vs = 0
    right_vs = 0
    left_max_level = -1
    right_max_level = -1
    max_vs_num = 0

    for tree in ka.kdata.field.trees.data
        for ps_data in tree
            hasproperty(ps_data, :vs_data) || continue
            vs_data = ps_data.vs_data
            vmax = maximum(vs_data.level)
            max_vs_num = max(max_vs_num, vs_data.vs_num)
            if ps_data.midpoint[1] < SPLIT_X
                left_ps += 1
                left_vs += vs_data.vs_num
                left_max_level = max(left_max_level, vmax)
            else
                right_ps += 1
                right_vs += vs_data.vs_num
                right_max_level = max(right_max_level, vmax)
            end
        end
    end

    return (; left_ps, right_ps, left_vs, right_vs,
            left_max_level, right_max_level, max_vs_num)
end

function global_vs_stats(ka)
    s = local_vs_stats(ka)
    left_ps = MPI.Allreduce(s.left_ps, MPI.SUM, MPI.COMM_WORLD)
    right_ps = MPI.Allreduce(s.right_ps, MPI.SUM, MPI.COMM_WORLD)
    left_vs = MPI.Allreduce(s.left_vs, MPI.SUM, MPI.COMM_WORLD)
    right_vs = MPI.Allreduce(s.right_vs, MPI.SUM, MPI.COMM_WORLD)
    left_max_level = MPI.Allreduce(s.left_max_level, MPI.MAX, MPI.COMM_WORLD)
    right_max_level = MPI.Allreduce(s.right_max_level, MPI.MAX, MPI.COMM_WORLD)
    max_vs_num = MPI.Allreduce(s.max_vs_num, MPI.MAX, MPI.COMM_WORLD)
    return (; left_ps, right_ps, left_vs, right_vs,
            total_ps = left_ps + right_ps,
            total_vs = left_vs + right_vs,
            left_max_level, right_max_level, max_vs_num)
end

function print_stats(label, stage, ka)
    rank = MPI.Comm_rank(MPI.COMM_WORLD)
    s = global_vs_stats(ka)
    if rank == 0
        @printf("[%s/%s] ps=%d vs=%d max_vs_num=%d left(maxL=%d, vs=%d) right(maxL=%d, vs=%d)\n",
                label, stage, s.total_ps, s.total_vs, s.max_vs_num,
                s.left_max_level, s.left_vs, s.right_max_level, s.right_vs)
    end
    return s
end

function write_case_summary(case_dir, label, initial_stats, final_stats, ka)
    MPI.Comm_rank(MPI.COMM_WORLD) == 0 || return nothing
    mkpath(case_dir)
    open(joinpath(case_dir, "case_summary.txt"), "w") do io
        println(io, "label = ", label)
        println(io, "sim_time = ", ka.kinfo.status.sim_time)
        println(io, "step = ", ka.kinfo.status.step)
        println(io, "mpi_size = ", MPI.Comm_size(MPI.COMM_WORLD))
        println(io, "geometry = ", ka.kinfo.config.geometry)
        println(io, "trees_num = ", ka.kinfo.config.trees_num)
        println(io, "vs_trees_num = ", ka.kinfo.config.vs_trees_num)
        println(io, "global_AMR_VS_MAXLEVEL = ",
                ka.kinfo.config.solver.AMR_VS_MAXLEVEL)
        println(io, "left_target_level = ", LEFT_VS_LEVEL)
        println(io, "right_target_level = ",
                occursin("L4_L3", label) ? RIGHT_VS_LEVEL : LEFT_VS_LEVEL)
        println(io, "initial_stats = ", initial_stats)
        println(io, "final_stats = ", final_stats)
    end
    return nothing
end

function run_case(label; heterogeneous = false)
    config = make_config()
    p4est, ka = initialize(config; prerefine_steps = 0)
    if heterogeneous
        changed_cells = cap_velocity_maxlevel!(ka)
        MPI.Comm_rank(MPI.COMM_WORLD) == 0 &&
            println("[$label/init] capped velocity grids in $changed_cells physical cells")
        amr_recover!(p4est, ka; topology_changed = true, velocity_changed = true)
    end

    initial_stats = print_stats(label, "init", ka)
    solve!(p4est, ka;
        animation = false,
        listen_for_save = false,
        status_check = false,
        progress = false,
        break_on_convergence = false,
    )
    final_stats = print_stats(label, "final", ka)

    case_dir = joinpath(OUTPUT_ROOT, label)
    write_case_summary(case_dir, label, initial_stats, final_stats, ka)
    save_result(p4est, ka; dir_path = case_dir)
    finalize!(p4est, ka)
    return case_dir
end

try
    rank = MPI.Comm_rank(MPI.COMM_WORLD)
    rank == 0 && println("Saving static VS mapping study under: ", OUTPUT_ROOT)
    if rank == 0
        println("Physical trees_num = ", PHYSICAL_TREES_NUM,
                ", MPI ranks = ", MPI.Comm_size(MPI.COMM_WORLD),
                ", max_sim_time = ", MAX_SIM_TIME)
    end

    grid_label = "Nx$(PHYSICAL_TREES_NUM[1])_Ny$(PHYSICAL_TREES_NUM[2])"
    l4_l4_dir = run_case("rp1d_static_L4_L4_$(grid_label)_t02"; heterogeneous = false)
    l4_l3_dir = run_case("rp1d_static_L4_L3_$(grid_label)_t02"; heterogeneous = true)

    if rank == 0
        println("Saved reference case: ", l4_l4_dir)
        println("Saved heterogeneous case: ", l4_l3_dir)
    end
finally
    MPI.Finalize()
end
