"""
To execute the code:
- from the main directory: mpiexecjl -n N --project=. julia example/cavity_sim.jl
- from the file directory: mpiexecjl -n N --project=.. julia cavity_sim.jl
"""

using KitAMR, MPI, CairoMakie

MPI.Init()
ps4est, DVM_data = KitAMR.init()
if MPI.Comm_rank(MPI.COMM_WORLD) == 0
    KitAMR.save_set(DVM_data)
end
KitAMR.save_result(DVM_data)
while true
    if DVM_data.global_data.adapt_step == 10
        KitAMR.update_slope!(DVM_data)
        KitAMR.update_gradmax!(DVM_data)
        KitAMR.PS_refine!(ps4est)
        KitAMR.PS_coarsen!(ps4est)
        KitAMR.PS_balance!(ps4est)
        if DVM_data.global_data.partition_step == 30
            KitAMR.PS_partition!(ps4est, DVM_data)
            DVM_data.global_data.partition_step = 0
        end
        if DVM_data.global_data.vs_adapt_step == 50
            KitAMR.vs_refine!(DVM_data)
            KitAMR.vs_coarsen!(DVM_data)
            DVM_data.global_data.vs_adapt_step = 0
        end
        KitAMR.update_ghost!(ps4est, DVM_data)
        KitAMR.update_neighbor!(ps4est, DVM_data)
        KitAMR.update_faces!(ps4est, DVM_data)
        DVM_data.global_data.adapt_step = 0
    end
    KitAMR.update_Δt!(DVM_data) # 798.125 μs (41515 allocations: 806.91 KiB)
    KitAMR.update_slope!(DVM_data) # 1.304 s (22645678 allocations: 1.43 GiB) -> final: 336.005 ms (32818 allocations: 1.35 MiB)
    KitAMR.DVM_data_exchange!(ps4est, DVM_data, Val(2)) # 0.121942
    KitAMR.update_flux!(DVM_data) # 4.437 s (174354747 allocations: 5.30 GiB)->3.268 s (75782427 allocations: 4.84 GiB) -> final: 540.574 ms (1301847 allocations: 1.42 GiB)->411.685 ms (1262009 allocations: 1.51 GiB)
    KitAMR.update_volume!(DVM_data) # 1.262 s (59874178 allocations: 1.55 GiB) ->469.264 ms (9070978 allocations: 623.28 MiB)-> final 50.184 ms (248053 allocations: 158.25 MiB)
    KitAMR.DVM_data_exchange!(ps4est, DVM_data, Val(1))
    DVM_data.global_data.adapt_step += 1
    DVM_data.global_data.vs_adapt_step += 1
    DVM_data.global_data.partition_step += 1
    KitAMR.save_result(DVM_data)
    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        @show DVM_data.global_data.gas.sim_time
        pp = PointerWrapper(ps4est)
        @show pp.local_num_quadrants[]
    end
    DVM_data.global_data.gas.sim_time > DVM_data.data_set.end_time && break
    MPI.Barrier(MPI.COMM_WORLD)
    if MPI.Comm_rank(MPI.COMM_WORLD) == MPI.Comm_size(MPI.COMM_WORLD) - 1
        ps_data = DVM_data.trees.data[end][end]
        @show ps_data.vs_data.vs_num
    end
end
KitAMR.save_VS_final(DVM_data)
meshes = KitAMR.collect_mesh(DVM_data)
solutions = KitAMR.collect_solution(ps4est, DVM_data)
# write_result!(ps4est,DVM_data,"write_test")
if MPI.Comm_rank(MPI.COMM_WORLD) == 0
    x, y = KitAMR.reshape_mesh(meshes)
    f = Figure()
    ax = Axis(f[1, 1])
    KitAMR.mesh_plot(x, y, ax)
    save("mesh.png", f)

    #=x, y, variable = KitAMR.reshape_solutions(solutions, DVM_data.global_data, :prim, 4)
    variable = 1 ./ variable
    f = Figure()
    ax = Axis(f[1, 1])
    co = contourf!(x, y, variable)
    Colorbar(f[1, 2], co)
    save("test_vs_adaptive.png", f)=#
end

KitAMR.finalize!(ps4est, DVM_data)

phase_cells = 0
for i in eachindex(DVM_data.trees.data)
    for j in eachindex(DVM_data.trees.data[i])
        global phase_cells
        phase_cells += DVM_data.trees.data[i][j].vs_data.vs_num
    end
end
phase_cells = MPI.Allreduce(phase_cells, +, MPI.COMM_WORLD)
if MPI.Comm_rank(MPI.COMM_WORLD) == 0
    @show phase_cells
end
