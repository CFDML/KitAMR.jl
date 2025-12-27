using KitAMR,MPI
include("./sphere_udf.jl")
MPI.Init()
config = KitAMR.read_config("./example/sphere/configure_sphere.txt")
ps4est,amr = KitAMR.init(config);
KitAMR.save_result(ps4est,amr)
KitAMR.listen_for_save!()
max_sim_time = 20.
nt = max_sim_time/amr.global_data.status.Δt+1.0 |> floor |> Int
for i in 1:nt
    if MPI.Comm_rank(MPI.COMM_WORLD)==0
        @show i
    end
    KitAMR.adaptive!(ps4est,amr;vs_interval=10,partition_interval=40)
    KitAMR.update_slope!(amr)
    KitAMR.slope_exchange!(ps4est, amr) 
    KitAMR.update_solid_cell!(amr)
    KitAMR.solid_exchange!(ps4est, amr)
    KitAMR.update_solid_neighbor!(amr)
    KitAMR.flux!(amr) 
    KitAMR.iterate!(amr) 
    KitAMR.data_exchange!(ps4est, amr)
    KitAMR.check_for_convergence(amr)&&break
    KitAMR.check!(i,ps4est,amr)
end
KitAMR.save_result(ps4est,amr)
KitAMR.finalize!(ps4est,amr)
MPI.Finalize()

