using KitAMR, MPI

MPI.Init()
ps4est, DVM_data = KitAMR.init();
if MPI.Comm_rank(MPI.COMM_WORLD) == 0
    KitAMR.save_set(DVM_data)
end

# minimal solution algorithm
begin
    KitAMR.update_slope!(DVM_data)
    KitAMR.update_gradmax!(DVM_data)
    KitAMR.PS_refine!(ps4est)
    KitAMR.PS_coarsen!(ps4est)
    KitAMR.PS_balance!(ps4est)
    KitAMR.PS_partition!(ps4est, DVM_data)
    KitAMR.update_ghost!(ps4est, DVM_data)
    KitAMR.update_neighbor!(ps4est, DVM_data)
    KitAMR.update_faces!(ps4est, DVM_data)
    KitAMR.update_Δt!(DVM_data)
    KitAMR.update_slope!(DVM_data)
    KitAMR.DVM_data_exchange!(ps4est, DVM_data, Val(2))
    KitAMR.update_flux!(DVM_data)
    KitAMR.update_volume!(DVM_data)
    KitAMR.DVM_data_exchange!(ps4est, DVM_data, Val(1))
    KitAMR.save_result(DVM_data)
    MPI.Barrier(MPI.COMM_WORLD)
end

KitAMR.save_VS_final(DVM_data)
meshes = KitAMR.collect_mesh(DVM_data)
solutions = KitAMR.collect_solution(ps4est, DVM_data)

# skip plotting
if MPI.Comm_rank(MPI.COMM_WORLD) == 0
    x, y = KitAMR.reshape_mesh(meshes)
    x, y, variable = KitAMR.reshape_solutions(solutions, DVM_data.global_data, :prim, 4)
end

KitAMR.finalize!(ps4est, DVM_data)
