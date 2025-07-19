using KitAMR,MPI
include("./cylinder_udf.jl")
MPI.Init()
config = KitAMR.read_config("./example/cylinder/configure_cylinder_Kn0p01.txt")
ps4est,amr = KitAMR.init(config);
# Partition for IB
KitAMR.ps_partition!(ps4est, amr)
KitAMR.update_ghost!(ps4est, amr)
KitAMR.update_neighbor!(ps4est, amr)
KitAMR.update_solid!(amr)
KitAMR.update_faces!(ps4est, amr)
KitAMR.listen_for_save!()
max_sim_time = 20.
nt = max_sim_time/amr.global_data.status.Δt+1.0 |> floor |> Int
for i in 1:nt
    KitAMR.adaptive!(ps4est,amr;partition_interval=200,ps_interval=100,vs_interval=200)
    KitAMR.update_slope!(amr;buffer_steps=100,i)
    KitAMR.slope_exchange!(ps4est, amr) 
    KitAMR.update_solid_cell!(amr)
    KitAMR.data_exchange!(ps4est, amr)
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

