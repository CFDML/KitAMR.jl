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
    AMR_PS_MAXLEVEL = 4,
    AMR_VS_MAXLEVEL = 5,
    AMR_PS_DYNAMIC = true,
    AMR_VS_DYNAMIC = true,
    AMR_VS_MODE = :haar,
    flux = CAIDVM,
    time_marching = CIP_Marching,
    max_sim_time = 0.25,
)

gas = Gas(;
    K = 1.0,
    Kn = 1e3,
    ω = 0.81,
    ωᵣ = 0.81,
)

output = Output(
    solver;
    vtk_celltype = [Pixel, Triangle],
    vs_vtk_celltype = Pixel,
    anim_dt = 0.01
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
        Domain(Period, 1), Domain(Period, 2),
        Domain(Period, 3), Domain(Period, 4),
    ],
    output = output,
    gas = gas,
    user_defined = udf,
)

p4est, ka = initialize(config)

solve!(p4est, ka)
save_result(p4est, ka)
finalize!(p4est, ka)
MPI.Finalize()
