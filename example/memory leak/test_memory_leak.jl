using KitAMR,MPI
include("../cylinder/cylinder_udf.jl")
MPI.Init()
config = KitAMR.read_config("./example/cylinder/configure_cylinder_Kn0p01.txt")
ps4est,amr = KitAMR.init(config);
max_sim_time = 20.
nt = max_sim_time/amr.global_data.status.Δt+1.0 |> floor |> Int
for i in 1:nt
    N = 0
    if MPI.Comm_rank(MPI.COMM_WORLD)==134
        i%100==0&&(@show i)
        sleep(0.1)
    end
    # KitAMR.update_ghost!(ps4est,amr)
    if i>1
        MPI.Barrier(MPI.COMM_WORLD)
    end
    KitAMR.update_neighbor!(ps4est, amr)
    # KitAMR.update_solid!(amr)
    # KitAMR.update_faces!(ps4est,amr)
    
end
KitAMR.finalize!(ps4est,amr)
MPI.Finalize()