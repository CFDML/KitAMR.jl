using KitAMR,MPI

MPI.Init()
config = KitAMR.read_config("configure_2D.txt")
ps4est,amr = KitAMR.init(config);

for i = 1:50
    if amr.global_data.status.ps_adapt_step == 10
        KitAMR.update_slope!(amr)
        KitAMR.update_gradmax!(amr)
        KitAMR.ps_refine!(ps4est,amr)
        KitAMR.ps_coarsen!(ps4est)
        KitAMR.ps_balance!(ps4est)
        if amr.global_data.status.partition_step == 20
            KitAMR.IB_quadid_update!(amr)
            KitAMR.ps_partition!(ps4est, amr)
            KitAMR.IB_update!(ps4est,amr)
            amr.global_data.status.partition_step = 0
        end
        if amr.global_data.status.vs_adapt_step == 40
            KitAMR.vs_refine!(amr)
            KitAMR.vs_coarsen!(amr)
            amr.global_data.status.partition_step!=0&&KitAMR.IB_structure_update!(amr)
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
    KitAMR.IB_data_exchange!(amr)
    KitAMR.update_solid_cell!(amr)
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

# IB_cells = amr.field.IB_cells[1]
# for i in eachindex(IB_cells.IB_nodes)
#     for j in eachindex(IB_cells.IB_nodes[i])
#         isa(IB_cells.IB_nodes[i][j],KitAMR.PS_Data)&&(IB_cells.IB_nodes[i][j].prim[end] = 0.5)
#     end
# end

function fieldvalues_fn(ps_data)
    isa(ps_data,KitAMR.InsideSolidData) && return [0.]
    ps_data.bound_enc<0 && return [-1.]
    return [1/ps_data.prim[end]]
end

KitAMR.write_VTK(ps4est,"testIB",["T"],fieldvalues_fn)
# KitAMR.finalize!(ps4est,amr)
MPI.Finalize()
# amr.field.aux_points[1]
# f = KitAMR.Figure()
# ax = KitAMR.Axis(f[1, 1])
# for i in eachindex(amr.field.aux_points[1])
#     KitAMR.scatter!(amr.field.aux_points[1][i]...)
# end
# KitAMR.save("test.png",f)

length(amr.field.boundary.solid_cells[1].ps_datas)
encs = [x.bound_enc for x in amr.field.boundary.solid_cells[1].ps_datas]
size(amr.field.trees.data[1][1].vs_data.df)
size(amr.field.boundary.solid_cells[1].ps_datas[1].vs_data.df)