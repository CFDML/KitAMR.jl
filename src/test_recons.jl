using KitAMR,MPI
MPI.Init()
config = KitAMR.read_config("configure_UGKS.txt")
ps4est,amr = KitAMR.init(config);
KitAMR.listen_for_save!()
for i in 1:100000
    if amr.global_data.status.ps_adapt_step == 10
        KitAMR.update_slope!(amr)
        KitAMR.update_gradmax!(amr)
        KitAMR.ps_refine!(ps4est,amr)
        KitAMR.ps_coarsen!(ps4est)
        KitAMR.ps_balance!(ps4est)
        if amr.global_data.status.vs_adapt_step == 10
            KitAMR.vs_refine!(amr)
            KitAMR.vs_coarsen!(amr)
            # KitAMR.update_solid_vsnum(amr)
            amr.global_data.status.vs_adapt_step=0
        end
        if amr.global_data.status.partition_step%40==0&&KitAMR.partition_check(ps4est)
            KitAMR.IB_quadid_update!(ps4est,amr)
            KitAMR.ps_partition!(ps4est, amr)
            KitAMR.IB_partition!(ps4est,amr)
            amr.global_data.status.partition_step = 0
        end
        amr.global_data.status.vs_adapt_step==0&&amr.global_data.status.partition_step!=0&&KitAMR.IB_structure_update!(amr)
        KitAMR.update_ghost!(ps4est, amr)
        KitAMR.update_neighbor!(ps4est, amr)
        KitAMR.update_faces!(ps4est, amr)
        amr.global_data.status.ps_adapt_step = 0
    end
    KitAMR.update_Δt!(amr)
    KitAMR.update_slope!(amr)
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
    KitAMR.check_for_save!(ps4est,amr)
end

KitAMR.save_result(ps4est,amr)
KitAMR.finalize!(ps4est,amr)
MPI.Finalize()