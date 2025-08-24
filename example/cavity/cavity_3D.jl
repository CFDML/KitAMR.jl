using KitAMR,MPI
MPI.Init()
config = KitAMR.read_config("./example/cavity/configure_cavity_3D.txt")
ps4est,amr = KitAMR.init(config);
KitAMR.listen_for_save!()
max_sim_time = 20.
nt = max_sim_time/amr.global_data.status.Δt+1.0 |> floor |> Int
for i in 1:nt
    if false
        @show i
        @time "adaptive" KitAMR.adaptive!(ps4est,amr;vs_interval=40,partition_interval=160)
        @time "slope" KitAMR.update_slope!(amr)
        @time "slope exchange" KitAMR.slope_exchange!(ps4est, amr) 
        # @time "slope Barrier" MPI.Barrier(MPI.COMM_WORLD)
        @time "solid cell" KitAMR.update_solid_cell!(amr)
        @time "solid exchange" KitAMR.solid_exchange!(ps4est, amr)
        # @time "solid Barrier" MPI.Barrier(MPI.COMM_WORLD)
        @time "solid neighbor" KitAMR.update_solid_neighbor!(amr)
        @time "flux" KitAMR.flux!(amr) 
        @time "iterate" KitAMR.iterate!(amr) # !Barrier! by Allreduce in Δt_comm()
        @time "data exchange" KitAMR.data_exchange!(ps4est, amr)
        KitAMR.check_for_convergence(amr)&&break
        @time "check" KitAMR.check!(i,ps4est,amr)
    else
        KitAMR.adaptive!(ps4est,amr;vs_interval=40,partition_interval=160)
        KitAMR.update_slope!(amr)
        KitAMR.slope_exchange!(ps4est, amr) 
        # MPI.Barrier(MPI.COMM_WORLD)
        KitAMR.update_solid_cell!(amr)
        KitAMR.solid_exchange!(ps4est, amr)
        # MPI.Barrier(MPI.COMM_WORLD)
        KitAMR.update_solid_neighbor!(amr)
        KitAMR.flux!(amr) 
        KitAMR.iterate!(amr) 
        KitAMR.data_exchange!(ps4est, amr)
        KitAMR.check_for_convergence(amr)&&break
        KitAMR.check!(i,ps4est,amr)
    end
end
KitAMR.save_result(ps4est,amr)
KitAMR.finalize!(ps4est,amr)
MPI.Finalize()

