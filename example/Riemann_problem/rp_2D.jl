using KitAMR, MPI

MPI.Init()

lambda_from_pressure(rho, p) = 0.5 * rho / p

function state_from_rho_u_v_p(rho, u, v, p)
    return [rho, u, v, lambda_from_pressure(rho, p)]
end

function Riemann_2D_init(midpoint, kinfo)
    x = midpoint[1]
    y = midpoint[2]

    if x >= 0.0 && y >= 0.0
        return state_from_rho_u_v_p(0.5313, 0.0, 0.0, 0.4)
    elseif x < 0.0 && y >= 0.0
        return state_from_rho_u_v_p(1.0, 0.7276, 0.0, 1.0)
    elseif x < 0.0 && y < 0.0
        return state_from_rho_u_v_p(0.8, 0.0, 0.0, 1.0)
    else
        return state_from_rho_u_v_p(1.0, 0.0, 0.7276, 1.0)
    end
end

solver = Solver(;
    DIM = 2,
    NDF = 2,
    CFL = 0.4,
    AMR_PS_MAXLEVEL = 5,
    AMR_VS_MAXLEVEL = 3,
    PS_DYNAMIC_AMR = true,
    ADAPT_COEFFI_PS = 0.3,
    VS_DYNAMIC_AMR = true,
    flux = CAIDVM,
    # flux = DVM,
    time_marching = CIP_Marching,
    # time_marching = Euler,
    # time_marching = CAIDVM_Marching,
    max_sim_time = 0.25,
)

gas = Gas(;
    K = 1.0,
    Kn = 0.001,
    ω = 0.81,
    ωᵣ = 0.81,
)

output = Output(
    solver;
)

udf = UDF(;
)

config = Configure(solver;
    geometry = [-0.5, 0.5, -0.5, 0.5],
    trees_num = [16, 16],
    quadrature = [-5.0, 5.0, -5.0, 5.0],
    vs_trees_num = [8, 8],
    IC = PCoordFn(Riemann_2D_init),
    domain = [
        Domain(UniformOutflow, 1), Domain(UniformOutflow, 2),
        Domain(UniformOutflow, 3), Domain(UniformOutflow, 4),
    ],
    output = output,
    gas = gas,
    user_defined = udf,
)

p4est, ka = initialize(config)

for _ in 1:5
    KitAMR.ps_adaptive_mesh_refinement!(p4est, ka; recursive = false)
    reinitialize_initial_condition!(ka)
    KitAMR.amr_recover!(p4est, ka)
end

KitAMR.execute_check!(p4est, ka)
listen_for_save!()

check_for_animsave!(p4est, ka)   # initial animation frame (no-op unless output.anim_dt>0)
i = 0
while !reached_max_time(ka)
    global i += 1
    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        @show i
    end
    adaptive_mesh_refinement!(p4est, ka; ps_interval = 40, vs_interval = 40, partition_interval = 40, vs_balance = false)
    limit_Δt!(ka)   # shrink Δt to land exactly on the next animation frame / max_sim_time
    slope!(p4est, ka)
    flux!(p4est, ka)
    iterate!(p4est, ka)
    check_for_animsave!(p4est, ka)   # write a frame if this step landed on a frame time
    check!(p4est, ka)
    check_for_convergence(ka) && break
end

KitAMR.execute_check!(p4est, ka)
save_result(p4est, ka)
finalize!(p4est, ka)
MPI.Finalize()
