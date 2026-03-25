using KitAMR,MPI
MPI.Init()
config = KitAMR.read_config("./example/cavity/configure_cavity_3D_32.txt")
p4est,ka = KitAMR.initialize_KitAMR(config);
max_sim_time = 20.
nt = max_sim_time/ka.kinfo.status.Δt+1.0 |> floor |> Int
for i in 1:200
    KitAMR.update_slope!(ka)
    KitAMR.slope_exchange!(p4est, ka) 
    KitAMR.flux!(ka) 
    KitAMR.iterate!(ka) 
    KitAMR.data_exchange!(p4est, ka)
end
KitAMR.finalize!(p4est,ka)
MPI.Finalize()

