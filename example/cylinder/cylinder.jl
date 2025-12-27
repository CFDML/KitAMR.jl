using KitAMR,MPI
include("./cylinder_udf.jl")
MPI.Init()
config = KitAMR.read_config("./example/cylinder/configure_cylinder.txt")
ps4est,amr = KitAMR.init(config);
KitAMR.listen_for_save!()
max_sim_time = 20.
nt = max_sim_time/amr.global_data.status.Δt+1.0 |> floor |> Int
# KitAMR.check_for_animsave!(ps4est,amr)
for i in 1:nt
    KitAMR.adaptive!(ps4est,amr;partition_interval=160)
    KitAMR.update_slope!(amr)
    KitAMR.slope_exchange!(ps4est, amr) 
    KitAMR.update_solid_cell!(amr)
    KitAMR.data_exchange!(ps4est, amr)
    KitAMR.update_solid_neighbor!(amr)
    KitAMR.slope_exchange!(ps4est, amr) 
    KitAMR.flux!(amr) 
    KitAMR.iterate!(amr) 
    KitAMR.data_exchange!(ps4est, amr)
    KitAMR.check_for_convergence(amr)&&break
    KitAMR.check!(i,ps4est,amr)
    # KitAMR.check_for_animsave!(ps4est,amr)
    # if amr.global_data.status.sim_time>3.0
    #     break
    # end
end
KitAMR.save_result(ps4est,amr)
KitAMR.finalize!(ps4est,amr)
MPI.Finalize()

