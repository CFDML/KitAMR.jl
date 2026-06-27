using KitAMR
using MPI
using Printf

const MAX_RANKS = 64

env_int(name::AbstractString, default::Integer) =
    something(tryparse(Int, get(ENV, name, "")), Int(default))

env_float(name::AbstractString, default::Real) =
    something(tryparse(Float64, get(ENV, name, "")), Float64(default))

function env_int3(name::AbstractString, default::NTuple{3,Int})
    raw = strip(get(ENV, name, ""))
    isempty(raw) && return collect(default)
    parts = split(replace(raw, "," => " "))
    if length(parts) == 1
        return fill(parse(Int, only(parts)), 3)
    elseif length(parts) == 3
        return parse.(Int, parts)
    end
    error("$name must be one integer or three comma/space separated integers.")
end

lambda_from_pressure(rho, p) = 0.5 * rho / p

state_from_rho_u_v_w_p(rho, u, v, w, p) =
    [rho, u, v, w, lambda_from_pressure(rho, p)]

function riemann_3d_init(midpoint, kinfo)
    x, y, z = midpoint
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

function stamp_qf_with_level!(p4est, ka)
    # Use the p8est quadrant level directly. Inferring the level from ps_data.ds
    # is more fragile after partition/balance and can disagree with the VTK
    # geometry traversal.
    config = ka.kinfo.config
    counts = zeros(Int, config.solver.AMR_PS_MAXLEVEL + 1)
    p_counts = pointer_from_objref(counts)
    GC.@preserve counts KitAMR.AMR_volume_iterate(p4est; user_data = p_counts) do ip, data, dp
        counts = unsafe_pointer_to_objref(data)
        ps_data = unsafe_pointer_to_objref(pointer(dp.ps_data))
        hasproperty(ps_data, :qf) || return nothing
        level = Int(ip.quad.level[])
        fill!(ps_data.qf, 0.0)
        ps_data.qf[1] = level
        counts[level + 1] += 1
    end
    return counts
end

function no_vs_output(ps_data, ka)
    return false
end

function riemann_interface_static_refine_flag(maxlevel::Integer, base_level::Integer, width_cells::Real)
    return (midpoint, ds, kinfo, level) -> begin
        level < base_level && return true
        level >= maxlevel && return false

        # Refine only the narrow cell layers cut by the initial Riemann interfaces.
        # The default width is slightly above 0.5 cell so that interfaces lying on
        # root-cell boundaries are still captured on both sides.
        for d in 1:3
            abs(midpoint[d]) <= width_cells * ds[d] && return true
        end
        return false
    end
end

function main()
    rank = MPI.Comm_rank(MPI.COMM_WORLD)
    nranks = MPI.Comm_size(MPI.COMM_WORLD)
    nranks <= MAX_RANKS ||
        error("rp_3D_level_vtk.jl is intended for at most $MAX_RANKS MPI ranks; got $nranks.")

    trees_num = env_int3("RP3D_VTK_TREES", (16, 16, 16))
    ps_maxlevel = env_int("RP3D_VTK_PS_MAXLEVEL", 4)
    base_level = min(env_int("RP3D_VTK_BASE_LEVEL", 0), ps_maxlevel)
    interface_width_cells = env_float("RP3D_VTK_INTERFACE_WIDTH_CELLS", 0.55)
    dynamic_steps = env_int("RP3D_VTK_DYNAMIC_STEPS", 0)
    vs_trees_num = env_int3("RP3D_VTK_VS_TREES", (1, 1, 1))
    v_bound = env_float("RP3D_VTK_V_BOUND", 1.0)
    outdir = get(ENV, "RP3D_VTK_OUTDIR", joinpath(@__DIR__, "rp_3D_level_vtk"))

    solver = Solver(;
        DIM = 3,
        NDF = 1,
        CFL = 0.4,
        AMR_PS_MAXLEVEL = ps_maxlevel,
        AMR_PS_DYNAMIC_MAXLEVEL = ps_maxlevel,
        AMR_VS_MAXLEVEL = 0,
        AMR_PS_DYNAMIC = false,
        AMR_VS_DYNAMIC = false,
        AMR_PS_THRES = env_float("RP3D_VTK_PS_THRES", 0.3),
        flux = CAIDVM,
        time_marching = CIP_Marching,
        max_sim_time = 0.0,
    )

    gas = Gas(;
        K = 1.0,
        Kn = 1e3,
        ω = 0.81,
        ωᵣ = 0.81,
    )

    output = Output(solver;
        vtk_celltype = Voxel,
        vs_vtk_celltype = Voxel,
        anim_dt = 0.0,
        vs_output_criterion = no_vs_output,
    )

    config = Configure(solver;
        geometry = [-0.5, 0.5, -0.5, 0.5, -0.5, 0.5],
        trees_num = trees_num,
        quadrature = [0.0, v_bound, 0.0, v_bound, 0.0, v_bound],
        vs_trees_num = vs_trees_num,
        IC = PCoordFn(riemann_3d_init),
        domain = [
            Domain(Period, 1), Domain(Period, 2),
            Domain(Period, 3), Domain(Period, 4),
            Domain(Period, 5), Domain(Period, 6),
        ],
        output = output,
        gas = gas,
        user_defined = UDF(;
            static_ps_refine_flag =
                riemann_interface_static_refine_flag(ps_maxlevel, base_level, interface_width_cells),
        ),
    )

    p4est = nothing
    ka = nothing
    try
        if rank == 0 && haskey(ENV, "RP3D_VTK_PREREFINE_STEPS")
            @printf("[rp3d_vtk] note: RP3D_VTK_PREREFINE_STEPS is ignored here; use RP3D_VTK_DYNAMIC_STEPS only if dynamic prerefinement is intentionally needed.\n")
        end
        rank == 0 && @printf(
            "[rp3d_vtk] ranks=%d trees=[%d,%d,%d] ps_lmax=%d base_level=%d interface_width_cells=%.3g dynamic_steps=%d vs_trees=[%d,%d,%d] outdir=%s\n",
            nranks, trees_num..., ps_maxlevel, base_level, interface_width_cells, dynamic_steps,
            vs_trees_num..., outdir)

        p4est, ka = initialize(config;
            prerefine_steps = dynamic_steps,
            prerefine_reinit_ic = true)

        local_counts = stamp_qf_with_level!(p4est, ka)
        global_counts = MPI.Allreduce(local_counts, MPI.SUM, MPI.COMM_WORLD)

        vtk_dir = joinpath(outdir, "vtk")
        rank == 0 && mkpath(vtk_dir)
        MPI.Barrier(MPI.COMM_WORLD)
        KitAMR.save_pvtu(joinpath(vtk_dir, "field"), p4est, ka, Voxel)
        MPI.Barrier(MPI.COMM_WORLD)

        if rank == 0
            @printf("[rp3d_vtk] total_physical_cells=%d\n", sum(global_counts))
            for level in 0:ps_maxlevel
                @printf("[rp3d_vtk] level=%d cells=%d\n", level, global_counts[level + 1])
            end
            @printf("[rp3d_vtk] wrote %s\n", joinpath(vtk_dir, "field.pvtu"))
            @printf("[rp3d_vtk] VTK cell data qf[1] stores the physical-space refinement level.\n")
        end
    finally
        if p4est !== nothing && ka !== nothing
            finalize!(p4est, ka)
        end
    end
end

MPI.Init()
try
    main()
finally
    MPI.Finalize()
end
