using KitAMR,MPI
MPI.Init()
config = KitAMR.read_config("./example/convergence_ps/configure_convergence_p8.txt")
ps4est,amr = KitAMR.init(config);
KitAMR.listen_for_save!()
max_sim_time = 20.
nt = max_sim_time/amr.global_data.status.Δt+1.0 |> floor |> Int
for i in 1:nt
    KitAMR.adaptive!(ps4est,amr;partition_interval=400,vs_interval=200)
    KitAMR.update_slope!(amr)
    KitAMR.slope_exchange!(ps4est, amr) 
    KitAMR.update_solid_cell!(amr)
    KitAMR.data_exchange!(ps4est, amr)
    KitAMR.update_solid_neighbor!(amr)
    KitAMR.flux!(amr) 
    KitAMR.iterate!(amr) 
    KitAMR.data_exchange!(ps4est, amr)
    KitAMR.check_for_convergence(amr)&&break
    KitAMR.check!(i,ps4est,amr)
    isnan(maximum(amr.global_data.status.residual.residual))&&break
end
KitAMR.save_result(ps4est,amr)
KitAMR.finalize!(ps4est,amr)
if MPI.Comm_rank(MPI.COMM_WORLD)==0
    run(`mv p ./example/convergence_ps/p/p8`)
end
MPI.Barrier(MPI.COMM_WORLD)

MPI.Finalize()

