using KitAMR,MPI
MPI.Init()
config = KitAMR.read_config("configure_convergence_64_v20.txt")
ps4est,amr = KitAMR.init(config);
KitAMR.listen_for_save!()
max_sim_time = 20.
nt = max_sim_time/amr.global_data.status.Δt+1.0 |> floor |> Int
for i in 1:nt
    KitAMR.adaptive!(ps4est,amr;partition_interval=160)
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
    run(`mv p ./test_convergence_vs/p20`)
end
MPI.Barrier(MPI.COMM_WORLD)


config = KitAMR.read_config("configure_convergence_64_v40.txt")
ps4est,amr = KitAMR.init(config);
max_sim_time = 20.
nt = max_sim_time/amr.global_data.status.Δt+1.0 |> floor |> Int
for i in 1:nt
    KitAMR.adaptive!(ps4est,amr;partition_interval=160)
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
    run(`mv p ./test_convergence_vs/p40`)
end
MPI.Barrier(MPI.COMM_WORLD)

config = KitAMR.read_config("configure_convergence_64_v80.txt")
ps4est,amr = KitAMR.init(config);
max_sim_time = 20.
nt = max_sim_time/amr.global_data.status.Δt+1.0 |> floor |> Int
for i in 1:nt
    KitAMR.adaptive!(ps4est,amr;partition_interval=160)
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
    run(`mv p ./test_convergence_vs/p80`)
end

config = KitAMR.read_config("configure_convergence_64_v200.txt")
ps4est,amr = KitAMR.init(config);
max_sim_time = 20.
nt = max_sim_time/amr.global_data.status.Δt+1.0 |> floor |> Int
for i in 1:nt
    KitAMR.adaptive!(ps4est,amr;partition_interval=160)
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
    run(`mv p ./test_convergence_vs/p200`)
end
MPI.Finalize()

