using KitAMR,MPI

MPI.Init()
config = KitAMR.read_config("configure_2D.txt")
ps4est,amr = KitAMR.init(config);

function fieldvalues_fn(ps_data)
    return [1/ps_data.prim[end]]
end

KitAMR.write_VTK(ps4est,"testIB",["T"],fieldvalues_fn)