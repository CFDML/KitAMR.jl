using KitAMR,MPI

MPI.Init()
config = KitAMR.read_config("configure_2D.txt")
ps4est,amr = KitAMR.init(config);


IB_cells = amr.field.IB_cells[1]
for i in eachindex(IB_cells.IB_nodes)
    for j in eachindex(IB_cells.IB_nodes[i])
        isa(IB_cells.IB_nodes[i][j],KitAMR.PS_Data)&&(IB_cells.IB_nodes[i][j].prim[end] = 0.5)
    end
end

function fieldvalues_fn(ps_data)
    isa(ps_data,KitAMR.InsideSolidData) && return [0.]
    return [1/ps_data.prim[end]]
end

KitAMR.write_VTK(ps4est,"testIB",["T"],fieldvalues_fn)
KitAMR.finalize!(ps4est,amr)
MPI.Finalize()
# amr.field.aux_points[1]
# f = KitAMR.Figure()
# ax = KitAMR.Axis(f[1, 1])
# for i in eachindex(amr.field.aux_points[1])
#     KitAMR.scatter!(amr.field.aux_points[1][i]...)
# end
# KitAMR.save("test.png",f)
