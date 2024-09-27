using KitAMR,MPI

MPI.Init()
config = KitAMR.read_config("configure_2D.txt")
ps4est,amr = KitAMR.init(config);

function fieldvalues_fn(ps_data)
    isa(ps_data,KitAMR.InsideSolidData) && return [0.]
    return [1/ps_data.prim[end]]
end

KitAMR.write_VTK(ps4est,"testIB",["T"],fieldvalues_fn)