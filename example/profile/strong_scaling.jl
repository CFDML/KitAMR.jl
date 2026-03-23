Pkg.add(NVTX)
using KitAMR, MPI, NVTX
MPI.Init()

solver = Solver(;
    DIM = 3,
    NDF = 1,
    AMR_PS_MAXLEVEL = 0,
    AMR_VS_MAXLEVEL = 0,
    PS_DYNAMIC_AMR = false,
    VS_DYNAMIC_AMR = false,
    flux = CAIDVM,
    time_marching = CAIDVM_Marching,
)
gas = Gas(; K = 0.0, Kn = 0.075, ω = 0.81, ωᵣ = 0.81)
output = Output(solver)
udf = UDF()
config = Configure(
    solver;
    geometry = [-0.5, 0.5, -0.5, 0.5, -0.5, 0.5],
    trees_num = [32, 32, 32],
    quadrature = [-5.0, 5.0, -5.0, 5.0, -5.0, 5.0],
    vs_trees_num = [24, 24, 24],
    IC = Uniform([1.0, 0.0, 0.0, 0.0, 1.0]),
    domain = [
        Domain(Maxwellian, 1, [1.0, 0.0, 0.0, 0.0, 1.0]),
        Domain(Maxwellian, 2, [1.0, 1.0*sqrt(5/6), 0.0, 0.0, 1.0]),
        Domain(Period, 3),
        Domain(Period, 4),
        Domain(Period, 5),
        Domain(Period, 6),
    ],
    output = output,
    gas = gas,
    user_defined = udf,
)

ps4est, amr = initialize_KitAMR(config);
KitAMR.listen_for_save!()
max_sim_time = 20.0
nt = max_sim_time/amr.global_data.status.Δt+1.0 |> floor |> Int
for i = 1:10
    # adaptive_mesh_refinement!(ps4est,amr)
    if MPI.Comm_rank(MPI.COMM_WORLD)==0
        @show i
    end
    NVTX.@range "update_slope" begin
        update_slope!(amr)
    end
    MPI.Barrier(MPI.COMM_WORLD)
    NVTX.@range "slope_exchange" begin
        slope_exchange!(ps4est, amr)
    end
    NVTX.@range "update_solid_cell" begin
        update_solid_cell!(amr)
    end
    NVTX.@range "solid_exchange" begin
        solid_exchange!(ps4est, amr)
    end
    NVTX.@range "update_solid_neighbor" begin
        update_solid_neighbor!(amr)
    end
    NVTX.@range "flux" begin
        flux!(amr)
    end
    NVTX.@range "iterate" begin
        iterate!(amr)
    end
    MPI.Barrier(MPI.COMM_WORLD)
    NVTX.@range "data_exchange" begin
        data_exchange!(ps4est, amr)
    end
    check_for_convergence(amr)&&break
    check!(i, ps4est, amr)
end
# save_result(ps4est,amr)
finalize!(ps4est, amr)
MPI.Finalize()
