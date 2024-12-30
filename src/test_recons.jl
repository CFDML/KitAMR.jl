using KitAMR,MPI
MPI.Init()
config = KitAMR.read_config("configure.txt")
ps4est,amr = KitAMR.init(config);
for i in 1:500
    KitAMR.update_Δt!(amr);
    KitAMR.update_slope!(amr);
    KitAMR.slope_exchange!(ps4est, amr) 
    KitAMR.flux!(amr) 
    KitAMR.update_volume!(amr) 
    KitAMR.IB_data_exchange!(amr)
    KitAMR.update_solid_cell!(amr)
    KitAMR.data_exchange!(ps4est, amr)
    amr.global_data.status.ps_adapt_step += 1
    amr.global_data.status.vs_adapt_step += 1
    amr.global_data.status.partition_step += 1
    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
	@show i
        @show amr.global_data.status.sim_time
        pp = KitAMR.PointerWrapper(ps4est)
        @show pp.local_num_quadrants[]
    end
    MPI.Barrier(MPI.COMM_WORLD)
    if MPI.Comm_rank(MPI.COMM_WORLD) == MPI.Comm_size(MPI.COMM_WORLD) - 1
        ps_data = amr.field.trees.data[end][end]
        @show ps_data.vs_data.vs_num
    end
end
KitAMR.save_result(ps4est,amr)
KitAMR.finalize!(ps4est,amr)
MPI.Finalize()