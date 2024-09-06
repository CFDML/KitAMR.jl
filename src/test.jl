using KitAMR,MPI,CairoMakie

MPI.Init()
config = KitAMR.read_config("configure.txt")
ps4est,amr = KitAMR.init(config);
for i = 1:100
    if amr.global_data.status.ps_adapt_step == 10
        KitAMR.update_slope!(amr)
        KitAMR.update_gradmax!(amr)
        KitAMR.ps_refine!(ps4est,amr)
        KitAMR.ps_coarsen!(ps4est)
        KitAMR.ps_balance!(ps4est)
        if amr.global_data.status.partition_step == 20
            KitAMR.ps_partition!(ps4est, amr)
            amr.global_data.status.partition_step = 0
        end
        if amr.global_data.status.vs_adapt_step == 40
            KitAMR.vs_refine!(amr)
            KitAMR.vs_coarsen!(amr)
            amr.global_data.status.vs_adapt_step = 0
        end
        KitAMR.update_ghost!(ps4est, amr)
        KitAMR.update_neighbor!(ps4est, amr)
        KitAMR.update_faces!(ps4est, amr)
        amr.global_data.status.ps_adapt_step = 0
    end
    KitAMR.update_Δt!(amr) 
    KitAMR.update_slope!(amr) 
    KitAMR.slope_exchange!(ps4est, amr) 
    KitAMR.update_flux!(amr) 
    # 2.540001 seconds (30.35 M allocations: 6.902 GiB, 11.05% gc time) vcat()
    # 5.969858 seconds (45.08 M allocations: 9.719 GiB, 9.98% gc time, 0.05% compilation time) Vcat()
    KitAMR.update_volume!(amr) 
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


solutions = KitAMR.collect_results(ps4est,amr)
dir_path = "./result/PS_result/"
!isdir(dir_path) && mkpath(dir_path)
KitAMR.JLD2.save_object(
        dir_path * "/" * string(MPI.Comm_rank(MPI.COMM_WORLD)) * ".jld",
        solutions
    )
if MPI.Comm_rank(MPI.COMM_WORLD)==0
    dir_path = "./result/SolverSet/"
    !isdir(dir_path) && mkpath(dir_path)
    KitAMR.JLD2.save_object(dir_path*"/SolverSet.jld",KitAMR.SolverSet(MPI.Comm_size(MPI.COMM_WORLD),amr.global_data))
end

# write_result!(ps4est,AMR_2D,"write_test")

 if MPI.Comm_rank(MPI.COMM_WORLD) == 0
     x, y, z, variable = KitAMR.reshape_solutions(solutions, amr.global_data, :prim, 5)
     vxz = variable[:,10,:]
     temp = 1 ./vxz
     f = Figure()
     ax = Axis(f[1, 1])
     co = contourf!(x, z, temp)
     Colorbar(f[1, 2], co)
     save("test.png", f)
 end

# a = rand(2000,2)
# b = rand(2000,2)
# c1 = @views CatView(a[1:1000,:],b[1:1000,:])
# c2 = @views vcat(a[1:1000,:],b[1:1000,:])
# @btime sum($c1)
# @btime sum($c2)

# using LazyArrays
# a = rand(2000,2)
# b = rand(2000,2)
# c1 = @views vcat(a[1:1000,:],b[1:1000,:])
# c2 = @views Vcat(a[1:1000,:],b[1:1000,:])
# d1 = @views c1[:,1]
# d2 = @views c2[:,1]
# e1 = @views c1[:,2]
# e2 = @views c2[:,2]
# @btime sum($d1)
# @btime sum($d2)
# @btime f1 = d1.+e1
# @btime f2 = d2.+e2
# f1 = d1.+e1
# f2 = d2.+e2
# t = Vector{Float64}
# @btime Vector{Float64}(f2)
# @btime 3.0.*$f1
# @btime 3.0.*$f2
