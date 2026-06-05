using KitAMR, MPI
include("./cylinder_udf.jl")
MPI.Init()

solver = Solver(;
    DIM = 2,
    NDF = 2,
    AMR_PS_MAXLEVEL = 7,
    AMR_VS_MAXLEVEL = 0,
    PS_DYNAMIC_AMR = false,
    VS_DYNAMIC_AMR = false,
    flux = CAIDVM,
    time_marching = CAIDVM_Marching,
)
gas = Gas(; K = 1.0, Kn = 0.1, ω = 0.81, ωᵣ = 0.81)
output = Output(solver;)
udf = UDF(; static_ps_refine_flag = shock_wave_region)
config = Configure(
    solver;
    geometry = [-16.0, 16.0, -16.0, 16.0],
    trees_num = [25, 25],
    quadrature = [-10.0, 10.0, -10.0, 10.0],
    vs_trees_num = [60, 60],
    IC = PCoordFn(cylinder_buffer_IC),
    domain = [
        Domain(SuperSonicInflow, 1, [1.0, 5.0*sqrt(5/6), 0.0, 1.0]),
        Domain(UniformOutflow, 2),
        Domain(UniformOutflow, 3),
        Domain(UniformOutflow, 4),
    ],
    IB = [Circle(Maxwellian, [0.0, 0.0], 1.0, true, 4.0, [1.0, 0.0, 0.0, 1.0])],
    output = output,
    gas = gas,
    user_defined = udf,
)


ps4est, amr = initialize_KitAMR(config);
listen_for_save!()
max_sim_time = 20.0
nt = max_sim_time/amr.global_data.status.Δt+1.0 |> floor |> Int
for i = 1:nt
    adaptive_mesh_refinement!(ps4est, amr; partition_interval = 160)
    update_slope!(amr)
    slope_exchange!(ps4est, amr)
    update_solid_cell!(amr)
    solid_exchange!(ps4est, amr)
    update_solid_neighbor!(amr)
    slope_exchange!(ps4est, amr)
    flux!(amr)
    iterate!(amr)
    data_exchange!(ps4est, amr)
    check_for_convergence(amr)&&break
    check!(i, ps4est, amr)
    # KitAMR.check_for_animsave!(ps4est,amr)
    # if amr.global_data.status.sim_time>3.0
    #     break
    # end
end
save_result(ps4est, amr)
finalize!(ps4est, amr)
MPI.Finalize()
