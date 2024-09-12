using KitAMR,MPI,CairoMakie

MPI.Init()
config = KitAMR.read_config("configure.txt")
ps4est,amr = KitAMR.init(config);


function fieldvalues_fn(ps_data)
    return [1/ps_data.prim[end]]
end

KitAMR.write_VTK(ps4est,"testT",["T"],fieldvalues_fn)

function fieldvalues_fn(vs_data::KitAMR.VS_Data{3})
    return [vs_data.df[:,1],vs_data.level]
end
KitAMR.write_vs_VTK(amr.field.trees.data[end][end].vs_data,amr,"testT_vs",["DistributionFunction,level"],fieldvalues_fn)