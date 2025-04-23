using KitAMR,MPI
include("./airfoil_udf.jl")
MPI.Init()
config = KitAMR.read_config("./example/airfoil/configure_airfoil.txt")
ps4est,amr = KitAMR.init(config);
KitAMR.listen_for_save!()
max_sim_time = 20.
nt = max_sim_time/amr.global_data.status.Δt+1.0 |> floor |> Int
# if MPI.Comm_rank(MPI.COMM_WORLD)==0
#     for i in 1:nt
#         KitAMR.adaptive!(ps4est,amr;partition_interval=160)
#         @time KitAMR.update_slope!(amr)
#         @time KitAMR.slope_exchange!(ps4est, amr) 
#         KitAMR.update_solid_cell!(amr)
#         @time KitAMR.data_exchange!(ps4est, amr)
#         KitAMR.update_solid_neighbor!(amr)
#         @time KitAMR.flux!(amr) 
#         @time KitAMR.iterate!(amr) 
#         @time KitAMR.data_exchange!(ps4est, amr)
#         KitAMR.check_for_convergence(amr)&&break
#         KitAMR.check!(i,ps4est,amr)
#         isnan(maximum(amr.global_data.status.residual.residual))&&break
#     end
# else
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
# end
KitAMR.save_result(ps4est,amr)
KitAMR.finalize!(ps4est,amr)
MPI.Finalize()

