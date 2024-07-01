include("DVM.jl")
MPI.Init()
ps4est, DVM_data = init();
MPI.Barrier(MPI.COMM_WORLD)
# save_result(DVM_data)
# while true
DVM_test_time = @elapsed begin
    for i = 1:200
        if DVM_data.global_data.adapt_step == 20
            update_slope!(DVM_data)
            update_gradmax!(DVM_data)
            PS_refine!(ps4est)
            PS_coarsen!(ps4est)
            PS_balance!(ps4est)
            if DVM_data.global_data.partition_step == 20
                PS_partition!(ps4est, DVM_data)
                DVM_data.global_data.partition_step = 0
            end
            if DVM_data.global_data.vs_adapt_step == 20
                vs_refine!(DVM_data)
                vs_coarsen!(DVM_data)
                DVM_data.global_data.vs_adapt_step = 0
            end
            update_ghost!(ps4est, DVM_data)
            update_neighbor!(ps4est, DVM_data)
            update_faces!(ps4est, DVM_data)
            DVM_data.global_data.adapt_step = 0
        end
        update_Δt!(DVM_data) # 798.125 μs (41515 allocations: 806.91 KiB)
        update_slope!(DVM_data) # 1.304 s (22645678 allocations: 1.43 GiB) -> final: 336.005 ms (32818 allocations: 1.35 MiB)
        DVM_data_exchange!(ps4est, DVM_data, Val(2)) # 0.121942
        update_flux!(DVM_data) # 4.437 s (174354747 allocations: 5.30 GiB)->3.268 s (75782427 allocations: 4.84 GiB) -> final: 540.574 ms (1301847 allocations: 1.42 GiB)->411.685 ms (1262009 allocations: 1.51 GiB)
        update_volume!(DVM_data) # 1.262 s (59874178 allocations: 1.55 GiB) ->469.264 ms (9070978 allocations: 623.28 MiB)-> final 50.184 ms (248053 allocations: 158.25 MiB)
        DVM_data_exchange!(ps4est, DVM_data, Val(1))
        DVM_data.global_data.adapt_step += 1
        DVM_data.global_data.vs_adapt_step += 1
        DVM_data.global_data.partition_step += 1
        # save_result(DVM_data)
        if MPI.Comm_rank(MPI.COMM_WORLD) == 0
            @show DVM_data.global_data.gas.sim_time
            pp = PointerWrapper(ps4est)
            @show pp.local_num_quadrants[]
            @show DVM_data.global_data.gradmax
            #@show data_size used_memory
        end
        # DVM_data.global_data.gas.sim_time>DVM_data.data_set.end_time && break
        MPI.Barrier(MPI.COMM_WORLD)
        if MPI.Comm_rank(MPI.COMM_WORLD) == MPI.Comm_size(MPI.COMM_WORLD) - 1
            ps_data = DVM_data.trees.data[end][end]
            @show ps_data.vs_data.vs_num
        end
    end
    MPI.Barrier(MPI.COMM_WORLD)
end
data_size = Base.summarysize(DVM_data)
data_size = MPI.Allreduce(data_size, +, MPI.COMM_WORLD)
data_size += p4est_memory_used(ps4est)
data_size /= 1024^3
GC.gc()
used_memory = Int(Sys.total_memory() - Sys.free_memory()) / 1024^3
if MPI.Comm_rank(MPI.COMM_WORLD) == 0
    @show data_size used_memory
end
# save_VS_final(DVM_data)
# save_set(DVM_data)

# meshes = collect_mesh(DVM_data)
# solutions = collect_solution(ps4est,DVM_data)
# # write_result!(ps4est,DVM_data,"write_test")
# if MPI.Comm_rank(MPI.COMM_WORLD) == 0
#     x,y = reshape_mesh(meshes)
#     f = Figure()
#     ax = Axis(f[1, 1])
#     mesh_plot(x, y,ax)
#     save("mesh.png", f)

#     x, y, variable = reshape_solutions(solutions, DVM_data.global_data, :prim, 4)
#     variable = 1 ./ variable
#     f = Figure()
#     ax = Axis(f[1, 1])
#     co = contourf!(x, y, variable)
#     Colorbar(f[1, 2], co)
#     save("test_vs_adaptive.png", f)
# end
finalize!(ps4est, DVM_data)
phase_cells = 0
for i in eachindex(DVM_data.trees.data)
    for j in eachindex(DVM_data.trees.data[i])
        global phase_cells
        phase_cells += DVM_data.trees.data[i][j].vs_data.vs_num
    end
end
phase_cells = MPI.Allreduce(phase_cells, +, MPI.COMM_WORLD)
if MPI.Comm_rank(MPI.COMM_WORLD) == 0
    @show DVM_test_time phase_cells
end
MPI.Finalize()
# end
# # mpiexecjl -n 6 --project=. julia ./main.jl  2592.10s user 308.70s system 552% cpu 8:45.11 total
