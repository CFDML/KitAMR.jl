using KitAMR,MPI
include("./X38_udf.jl")
MPI.Init()
config = KitAMR.read_config("./example/X38/configure_X38.txt")
ps4est,amr = KitAMR.init(config);
KitAMR.listen_for_save!()
max_sim_time = 20.
nt = max_sim_time/amr.global_data.status.Δt+1.0 |> floor |> Int
for i in 1:nt
    if MPI.Comm_rank(MPI.COMM_WORLD)==0
        @show i
    end
    KitAMR.adaptive!(ps4est,amr;ps_interval = 40, vs_interval=40,partition_interval=40)
    KitAMR.update_slope!(amr)
    KitAMR.slope_exchange!(ps4est, amr) 
    KitAMR.update_solid_cell!(amr)
    KitAMR.solid_exchange!(ps4est, amr)
    KitAMR.update_solid_neighbor!(amr)
    KitAMR.flux!(amr) 
    KitAMR.iterate!(amr) 
    KitAMR.data_exchange!(ps4est, amr)
    KitAMR.check_for_convergence(amr)&&break
    KitAMR.check!(i,ps4est,amr)
    # if i==400
    #     KitAMR.save_result(ps4est,amr)
    # end
    # if i==100
    #     trees = amr.field.trees.data
    #     for tree in trees
    #         for ps_data in tree
    #             isa(ps_data,KitAMR.InsideSolidData)&&continue
    #             point = [0.70654296875, 0.02734375, 0.18359375]
    #             if KitAMR.norm(point-ps_data.midpoint)<1e-12
    #                 KitAMR.write_vs_VTK(ps_data.vs_data.df,ps_data.vs_data,amr,"./X38_singular_vs_step"*string(i),["df"],KitAMR.fieldvalues_fn)
    #             end
    #             # point = [0.70654296875, 0.02734375, -0.18359375]
    #             # if KitAMR.norm(point-ps_data.midpoint)<1e-12
    #             #     KitAMR.write_vs_VTK(ps_data.vs_data.df,ps_data.vs_data,amr,"./X38_symetric_vs_step132",["df"],KitAMR.fieldvalues_fn)
    #             # end
    #         end
    #     end
    # end
end
# trees = amr.field.trees.data
# for tree in trees
#     for ps_data in tree
#         isa(ps_data,KitAMR.InsideSolidData)&&continue
#         point = [0.70654296875, 0.02734375, 0.18359375]
#         if KitAMR.norm(point-ps_data.midpoint)<1e-12
#             # KitAMR.write_vs_VTK(ps_data.vs_data.df,ps_data.vs_data,amr,"./X38_singular_vs_step132",["df"],KitAMR.fieldvalues_fn)
#             neighbors = ps_data.neighbor.data
#             for j in 1:6
#                 if isa(neighbors[j][1],KitAMR.SolidNeighbor)
#                     solid_neighbor = neighbors[j][1]
#                     @show solid_neighbor.aux_point solid_neighbor.normal solid_neighbor.midpoint
#                 end
#             end
#         end

#         point = [0.70654296875, 0.02734375, -0.18359375]
#         if KitAMR.norm(point-ps_data.midpoint)<1e-12
#             # KitAMR.write_vs_VTK(ps_data.vs_data.df,ps_data.vs_data,amr,"./X38_singular_vs_step132",["df"],KitAMR.fieldvalues_fn)
#             neighbors = ps_data.neighbor.data
#             for j in 1:6
#                 if isa(neighbors[j][1],KitAMR.SolidNeighbor)
#                     solid_neighbor = neighbors[j][1]
#                     @show solid_neighbor.aux_point solid_neighbor.normal solid_neighbor.midpoint
#                 end
#             end
#         end
#     end
# end
KitAMR.save_result(ps4est,amr)
KitAMR.finalize!(ps4est,amr)
MPI.Finalize()

