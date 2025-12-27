using KitAMR,MPI
include("./convergence_ps_udf.jl")
MPI.Init()


# config = KitAMR.read_config("./example/artificial_interpolation/configure_convergence_p16.txt")
# ps4est,amr = KitAMR.init(config);
# max_sim_time = 20.
# nt = max_sim_time/amr.global_data.status.Δt+1.0 |> floor |> Int
# for i in 1:1
#     KitAMR.adaptive!(ps4est,amr;partition_interval=400,vs_interval=200)
#     KitAMR.update_slope!(amr)
#     KitAMR.slope_exchange!(ps4est, amr) 
#     KitAMR.update_solid_cell!(amr)
#     KitAMR.data_exchange!(ps4est, amr)
# end
# KitAMR.save_result(ps4est,amr;dir_path="result_interp_p16")
# KitAMR.finalize!(ps4est,amr)
# MPI.Barrier(MPI.COMM_WORLD)

# config = KitAMR.read_config("./example/artificial_interpolation/configure_convergence_p32.txt")
# ps4est,amr = KitAMR.init(config);
# max_sim_time = 20.
# nt = max_sim_time/amr.global_data.status.Δt+1.0 |> floor |> Int
# for i in 1:1
#     KitAMR.adaptive!(ps4est,amr;partition_interval=400,vs_interval=200)
#     KitAMR.update_slope!(amr)
#     KitAMR.slope_exchange!(ps4est, amr) 
#     KitAMR.update_solid_cell!(amr)
#     KitAMR.data_exchange!(ps4est, amr)
# end
# KitAMR.save_result(ps4est,amr;dir_path="result_interp_p32")
# KitAMR.finalize!(ps4est,amr)
# MPI.Barrier(MPI.COMM_WORLD)

config = KitAMR.read_config("./example/artificial_interpolation/configure_convergence_p64.txt")
ps4est,amr = KitAMR.init(config);
max_sim_time = 20.
nt = max_sim_time/amr.global_data.status.Δt+1.0 |> floor |> Int
for i in 1:1
    KitAMR.adaptive!(ps4est,amr;partition_interval=400,vs_interval=200)
    KitAMR.update_slope!(amr)
    KitAMR.slope_exchange!(ps4est, amr) 
    KitAMR.update_solid_cell!(amr)
    KitAMR.data_exchange!(ps4est, amr)
end
KitAMR.save_result(ps4est,amr;dir_path="result_interp_p64")
KitAMR.finalize!(ps4est,amr)
MPI.Barrier(MPI.COMM_WORLD)

config = KitAMR.read_config("./example/artificial_interpolation/configure_convergence_p128.txt")
ps4est,amr = KitAMR.init(config);
max_sim_time = 20.
nt = max_sim_time/amr.global_data.status.Δt+1.0 |> floor |> Int
for i in 1:1
    KitAMR.adaptive!(ps4est,amr;partition_interval=400,vs_interval=200)
    KitAMR.update_slope!(amr)
    KitAMR.slope_exchange!(ps4est, amr) 
    KitAMR.update_solid_cell!(amr)
    KitAMR.data_exchange!(ps4est, amr)
end
KitAMR.save_result(ps4est,amr;dir_path="result_interp_p128")
KitAMR.finalize!(ps4est,amr)
MPI.Barrier(MPI.COMM_WORLD)

config = KitAMR.read_config("./example/artificial_interpolation/configure_convergence_p256.txt")
ps4est,amr = KitAMR.init(config);
max_sim_time = 20.
nt = max_sim_time/amr.global_data.status.Δt+1.0 |> floor |> Int
for i in 1:1
    KitAMR.adaptive!(ps4est,amr;partition_interval=400,vs_interval=200)
    KitAMR.update_slope!(amr)
    KitAMR.slope_exchange!(ps4est, amr) 
    KitAMR.update_solid_cell!(amr)
    KitAMR.data_exchange!(ps4est, amr)
end
KitAMR.save_result(ps4est,amr;dir_path="result_interp_p256")
KitAMR.finalize!(ps4est,amr)
MPI.Barrier(MPI.COMM_WORLD)

MPI.Finalize()

