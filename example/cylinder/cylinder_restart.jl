using KitAMR, MPI
include("./cylinder_udf.jl")
MPI.Init()
p4est,ka = restart("cylinder_restart")
ka.kinfo.status.residual.residual .= Inf
ka.kinfo.status.sim_time = 0.0
solve!(p4est, ka;vs_balance=true,ps_interval = 40, vs_interval = 40, partition_interval = 40)
save_result(p4est,ka)
save_for_restart(p4est,ka;dir_path="cylinder_restart")
finalize!(p4est,ka)
MPI.Finalize()