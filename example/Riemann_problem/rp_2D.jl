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

gas = Gas(; K = 1.0, Kn = 0.001, ω = 0.81, ωᵣ = 0.81)

output = Output(solver;)

udf = UDF(;)

config = Configure(
    solver;
    geometry = [-0.5, 0.5, -0.5, 0.5],
    trees_num = [16, 16],
    quadrature = [-5.0, 5.0, -5.0, 5.0],
    vs_trees_num = [8, 8],
    IC = PCoordFn(Riemann_2D_init),
    domain = [
        Domain(UniformOutflow, 1),
        Domain(UniformOutflow, 2),
        Domain(UniformOutflow, 3),
        Domain(UniformOutflow, 4),
    ],
    output = output,
    gas = gas,
    user_defined = udf,
)

p4est, ka = initialize(config; prerefine_steps = 5, prerefine_reinit_ic = true)

solve!(p4est, ka; ps_interval = 40, vs_interval = 40, partition_interval = 40)
save_result(p4est, ka)
finalize!(p4est, ka)
MPI.Finalize()
