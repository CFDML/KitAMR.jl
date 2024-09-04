include("DVM.jl")
MPI.Init()
ps4est, AMR_2D = init();
MPI.Barrier(MPI.COMM_WORLD)
# save_result(AMR_2D)
# while true
DVM_test_time = @elapsed begin
    for i = 1:200
        if AMR_2D.global_data.adapt_step == 20
            update_slope!(AMR_2D)
            update_gradmax!(AMR_2D)
            PS_refine!(ps4est)
            PS_coarsen!(ps4est)
            PS_balance!(ps4est)
            if AMR_2D.global_data.partition_step == 20
                PS_partition!(ps4est, AMR_2D)
                AMR_2D.global_data.partition_step = 0
            end
            if AMR_2D.global_data.vs_adapt_step == 20
                vs_refine!(AMR_2D)
                vs_coarsen!(AMR_2D)
                AMR_2D.global_data.vs_adapt_step = 0
            end
            update_ghost!(ps4est, AMR_2D)
            update_neighbor!(ps4est, AMR_2D)
            update_faces!(ps4est, AMR_2D)
            AMR_2D.global_data.adapt_step = 0
        end
        update_Δt!(AMR_2D) # 798.125 μs (41515 allocations: 806.91 KiB)
        update_slope!(AMR_2D) # 1.304 s (22645678 allocations: 1.43 GiB) -> final: 336.005 ms (32818 allocations: 1.35 MiB)
        AMR_exchange!(ps4est, AMR_2D, Val(2)) # 0.121942
        update_flux!(AMR_2D) # 4.437 s (174354747 allocations: 5.30 GiB)->3.268 s (75782427 allocations: 4.84 GiB) -> final: 540.574 ms (1301847 allocations: 1.42 GiB)->411.685 ms (1262009 allocations: 1.51 GiB)
        update_volume!(AMR_2D) # 1.262 s (59874178 allocations: 1.55 GiB) ->469.264 ms (9070978 allocations: 623.28 MiB)-> final 50.184 ms (248053 allocations: 158.25 MiB)
        AMR_exchange!(ps4est, AMR_2D, Val(1))
        AMR_2D.global_data.adapt_step += 1
        AMR_2D.global_data.vs_adapt_step += 1
        AMR_2D.global_data.partition_step += 1
        # save_result(AMR_2D)
        if MPI.Comm_rank(MPI.COMM_WORLD) == 0
            @show AMR_2D.global_data.gas.sim_time
            pp = PointerWrapper(ps4est)
            @show pp.local_num_quadrants[]
            @show AMR_2D.global_data.gradmax
            #@show data_size used_memory
        end
        # AMR_2D.global_data.gas.sim_time>AMR_2D.data_set.end_time && break
        MPI.Barrier(MPI.COMM_WORLD)
        if MPI.Comm_rank(MPI.COMM_WORLD) == MPI.Comm_size(MPI.COMM_WORLD) - 1
            ps_data = AMR_2D.trees.data[end][end]
            @show ps_data.vs_data.vs_num
        end
    end
    MPI.Barrier(MPI.COMM_WORLD)
end
data_size = Base.summarysize(AMR_2D)
data_size = MPI.Allreduce(data_size, +, MPI.COMM_WORLD)
data_size += p4est_memory_used(ps4est)
data_size /= 1024^3
GC.gc()
used_memory = Int(Sys.total_memory() - Sys.free_memory()) / 1024^3
if MPI.Comm_rank(MPI.COMM_WORLD) == 0
    @show data_size used_memory
end
# save_VS_final(AMR_2D)
# save_set(AMR_2D)

# meshes = collect_mesh(AMR_2D)
# solutions = collect_solution(ps4est,AMR_2D)
# # write_result!(ps4est,AMR_2D,"write_test")
# if MPI.Comm_rank(MPI.COMM_WORLD) == 0
#     x,y = reshape_mesh(meshes)
#     f = Figure()
#     ax = Axis(f[1, 1])
#     mesh_plot(x, y,ax)
#     save("mesh.png", f)

#     x, y, variable = reshape_solutions(solutions, AMR_2D.global_data, :prim, 4)
#     variable = 1 ./ variable
#     f = Figure()
#     ax = Axis(f[1, 1])
#     co = contourf!(x, y, variable)
#     Colorbar(f[1, 2], co)
#     save("test_vs_adaptive.png", f)
# end
finalize!(ps4est, AMR_2D)
phase_cells = 0
for i in eachindex(AMR_2D.trees.data)
    for j in eachindex(AMR_2D.trees.data[i])
        global phase_cells
        phase_cells += AMR_2D.trees.data[i][j].vs_data.vs_num
    end
end
phase_cells = MPI.Allreduce(phase_cells, +, MPI.COMM_WORLD)
if MPI.Comm_rank(MPI.COMM_WORLD) == 0
    @show DVM_test_time phase_cells
end
MPI.Finalize()
# end
# # mpiexecjl -n 6 --project=. julia ./main.jl  2592.10s user 308.70s system 552% cpu 8:45.11 total
