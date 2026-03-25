using KitAMR,MPI
include("./X38_udf.jl")
MPI.Init()
config = KitAMR.read_config("./example/X38/configure_X38.txt")
p4est,ka = KitAMR.initialize_KitAMR(config);
KitAMR.listen_for_save!()
max_sim_time = 20.
nt = max_sim_time/ka.kinfo.status.Δt+1.0 |> floor |> Int
for i in 1:nt
    if MPI.Comm_rank(MPI.COMM_WORLD)==0
        @show i
    end
    KitAMR.adaptive_mesh_refinement!(p4est,ka;ps_interval = 40, vs_interval=40,partition_interval=40)
    KitAMR.update_slope!(ka)
    KitAMR.slope_exchange!(p4est, ka) 
    KitAMR.update_solid_cell!(ka)
    KitAMR.solid_exchange!(p4est, ka)
    KitAMR.update_solid_neighbor!(ka)
    KitAMR.flux!(ka) 
    KitAMR.iterate!(ka) 
    KitAMR.data_exchange!(p4est, ka)
    KitAMR.check_for_convergence(ka)&&break
    KitAMR.check!(p4est,ka)
end
KitAMR.save_result(p4est,ka)
KitAMR.finalize!(p4est,ka)
MPI.Finalize()

