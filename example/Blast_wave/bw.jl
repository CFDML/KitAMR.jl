using KitAMR, MPI

MPI.Init()

function blast_wave_init(midpoint, kinfo)
    geometry = kinfo.config.geometry
    Lx = geometry[2] - geometry[1]
    Ly = geometry[4] - geometry[3]
    x0 = 0.5 * (geometry[1] + geometry[2])
    y0 = 0.5 * (geometry[3] + geometry[4])
    dx = abs(midpoint[1] - x0)
    dy = abs(midpoint[2] - y0)
    dx = min(dx, Lx - dx)
    dy = min(dy, Ly - dy)
    r = sqrt(dx^2 + dy^2)

    rho = 1.0
    ux = 0.0
    uy = 0.0

    T_ambient = 1.0
    T_peak = 80.0
    radius = 0.08
    width = 0.01
    hot_fraction = 0.5 * (1.0 - tanh((r - radius) / width))
    T = T_ambient + (T_peak - T_ambient) * hot_fraction

    return [rho, ux, uy, 1.0 / T]
end

solver = Solver(;
    DIM = 2,
    NDF = 2,
    CFL = 0.35,
    AMR_PS_MAXLEVEL = 4,
    AMR_VS_MAXLEVEL = 5,
    PS_DYNAMIC_AMR = true,
    ADAPT_COEFFI_PS = 0.8,
    VS_DYNAMIC_AMR = true,
    flux = CAIDVM,
    time_marching = CIP_Marching,
)

gas = Gas(;
    K = 1.0,
    Kn = 1e-3,
    ω = 0.81,
    ωᵣ = 0.81,
)

output = Output(
    solver;
)

udf = UDF(;
)

config = Configure(solver;
    geometry = [0.0, 1.0, 0.0, 1.0],
    trees_num = [16, 16],
    quadrature = [-30.0, 30.0, -30.0, 30.0],
    vs_trees_num = [8, 8],
    IC = PCoordFn(blast_wave_init),
    domain = [
        Domain(Period, 1), Domain(Period, 2),
        Domain(Period, 3), Domain(Period, 4),
    ],
    output = output,
    gas = gas,
    user_defined = udf,
)

p4est, ka = initialize(config)
slope!(p4est, ka)

for _ in 1:3
    KitAMR.ps_adaptive_mesh_refinement!(p4est, ka; recursive = false)
    reinitialize_initial_condition!(ka)
    KitAMR.amr_recover!(p4est, ka)
end

KitAMR.execute_check!(p4est, ka)
listen_for_save!()

max_sim_time = 0.25
nt = max_sim_time / ka.kinfo.status.Δt + 1.0 |> floor |> Int
for i in 1:nt
    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        @show i
    end
    adaptive_mesh_refinement!(p4est, ka; ps_interval = 20, vs_interval = 40, partition_interval = 40, vs_balance = false)
    slope!(p4est, ka)
    flux!(p4est, ka)
    iterate!(p4est, ka)
    check!(p4est, ka)
    check_for_convergence(ka) && break
end

KitAMR.execute_check!(p4est, ka)
save_result(p4est, ka)
finalize!(p4est, ka)
MPI.Finalize()
