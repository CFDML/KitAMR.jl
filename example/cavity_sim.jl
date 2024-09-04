"""
To execute the code (assuming 16 cores are going to be used):
- from the main directory: mpiexecjl -n 16 --project=. julia example/cavity_sim.jl
- from the file directory: mpiexecjl -n 16 --project=.. julia cavity_sim.jl
"""

using KitAMR, MPI, CairoMakie

MPI.Init()
ps4est, AMR_2D = KitAMR.init()
if MPI.Comm_rank(MPI.COMM_WORLD) == 0
    KitAMR.save_set(AMR_2D)
end
KitAMR.save_result(AMR_2D)
while true
    if AMR_2D.global_data.adapt_step == 10
        KitAMR.update_slope!(AMR_2D)
        KitAMR.update_gradmax!(AMR_2D)
        KitAMR.PS_refine!(ps4est)
        KitAMR.PS_coarsen!(ps4est)
        KitAMR.PS_balance!(ps4est)
        if AMR_2D.global_data.partition_step == 30
            KitAMR.PS_partition!(ps4est, AMR_2D)
            AMR_2D.global_data.partition_step = 0
        end
        if AMR_2D.global_data.vs_adapt_step == 50
            KitAMR.vs_refine!(AMR_2D)
            KitAMR.vs_coarsen!(AMR_2D)
            AMR_2D.global_data.vs_adapt_step = 0
        end
        KitAMR.update_ghost!(ps4est, AMR_2D)
        KitAMR.update_neighbor!(ps4est, AMR_2D)
        KitAMR.update_faces!(ps4est, AMR_2D)
        AMR_2D.global_data.adapt_step = 0
    end
    KitAMR.update_Δt!(AMR_2D) # 798.125 μs (41515 allocations: 806.91 KiB)
    KitAMR.update_slope!(AMR_2D) # 1.304 s (22645678 allocations: 1.43 GiB) -> final: 336.005 ms (32818 allocations: 1.35 MiB)
    KitAMR.AMR_exchange!(ps4est, AMR_2D, Val(2)) # 0.121942
    KitAMR.update_flux!(AMR_2D) # 4.437 s (174354747 allocations: 5.30 GiB)->3.268 s (75782427 allocations: 4.84 GiB) -> final: 540.574 ms (1301847 allocations: 1.42 GiB)->411.685 ms (1262009 allocations: 1.51 GiB)
    KitAMR.update_volume!(AMR_2D) # 1.262 s (59874178 allocations: 1.55 GiB) ->469.264 ms (9070978 allocations: 623.28 MiB)-> final 50.184 ms (248053 allocations: 158.25 MiB)
    KitAMR.AMR_exchange!(ps4est, AMR_2D, Val(1))
    AMR_2D.global_data.adapt_step += 1
    AMR_2D.global_data.vs_adapt_step += 1
    AMR_2D.global_data.partition_step += 1
    KitAMR.save_result(AMR_2D)
    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        @show AMR_2D.global_data.gas.sim_time
        pp = PointerWrapper(ps4est)
        @show pp.local_num_quadrants[]
    end
    AMR_2D.global_data.gas.sim_time > AMR_2D.data_set.end_time && break
    MPI.Barrier(MPI.COMM_WORLD)
    if MPI.Comm_rank(MPI.COMM_WORLD) == MPI.Comm_size(MPI.COMM_WORLD) - 1
        ps_data = AMR_2D.trees.data[end][end]
        @show ps_data.vs_data.vs_num
    end
end
KitAMR.save_VS_final(AMR_2D)
meshes = KitAMR.collect_mesh(AMR_2D)
solutions = KitAMR.collect_solution(ps4est, AMR_2D)

if MPI.Comm_rank(MPI.COMM_WORLD) == 0
    x, y = KitAMR.reshape_mesh(meshes)
    f = Figure()
    ax = Axis(f[1, 1])
    KitAMR.mesh_plot(x, y, ax)
    save("mesh.png", f)

    x, y, variable = KitAMR.reshape_solutions(solutions, AMR_2D.global_data, :prim, 4)
    variable = 1 ./ variable
    f = Figure()
    ax = Axis(f[1, 1])
    co = contourf!(x, y, variable)
    Colorbar(f[1, 2], co)
    save("test_vs_adaptive.png", f)
end

KitAMR.finalize!(ps4est, AMR_2D)

phase_cells = 0
for i in eachindex(AMR_2D.trees.data)
    for j in eachindex(AMR_2D.trees.data[i])
        global phase_cells
        phase_cells += AMR_2D.trees.data[i][j].vs_data.vs_num
    end
end
phase_cells = MPI.Allreduce(phase_cells, +, MPI.COMM_WORLD)
if MPI.Comm_rank(MPI.COMM_WORLD) == 0
    @show phase_cells
end
