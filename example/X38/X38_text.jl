using KitAMR, MPI
include("./X38_udf.jl")
MPI.Init()
config = KitAMR.read_config("./example/X38/configure_X38.txt")
p4est, ka = KitAMR.initialize(config);
KitAMR.listen_for_save!()
i = 0
while !KitAMR.reached_max_time(ka)   # max_sim_time is read from the config file
    global i += 1
    if MPI.Comm_rank(MPI.COMM_WORLD)==0
        @show i
    end
    KitAMR.adaptive_mesh_refinement!(
        p4est,
        ka;
        ps_interval = 40,
        vs_interval = 40,
        partition_interval = 40,
    )
    KitAMR.limit_Δt!(ka)   # shrink Δt to land exactly on the next animation frame / max_sim_time
    KitAMR.update_slope!(ka)
    KitAMR.slope_exchange!(p4est, ka)
    KitAMR.update_solid_cell!(ka)
    KitAMR.solid_exchange!(p4est, ka)
    KitAMR.update_solid_neighbor!(ka)
    KitAMR.flux!(ka)
    KitAMR.iterate!(ka)
    KitAMR.data_exchange!(p4est, ka)
    check_for_animsave!(p4est, ka)   # write a frame if this step landed on a frame time
    KitAMR.check_for_convergence(ka)&&break
    KitAMR.check!(p4est, ka)
end
KitAMR.save_result(p4est, ka)
KitAMR.finalize!(p4est, ka)
MPI.Finalize()
