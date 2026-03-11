using KitAMR,MPI
MPI.Init()
config = KitAMR.read_config("./example/cavity/configure_cavity_3D_32.txt")
ps4est,amr = KitAMR.initialize_KitAMR(config);
max_sim_time = 20.
nt = max_sim_time/amr.global_data.status.Δt+1.0 |> floor |> Int
for i in 1:200
    KitAMR.update_slope!(amr)
    KitAMR.slope_exchange!(ps4est, amr) 
    KitAMR.flux!(amr) 
    KitAMR.iterate!(amr) 
    KitAMR.data_exchange!(ps4est, amr)
end
KitAMR.finalize!(ps4est,amr)
MPI.Finalize()

