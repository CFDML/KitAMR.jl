using KitAMR,MPI

MPI.Init()
config = KitAMR.read_config("configure_Cylinder.txt")
ps4est,amr = KitAMR.init(config);

fp = PointerWrapper(ps4est)
cp = fp.connectivity
cp.num_corners[]
cp.num_trees[]
cp.num_vertices[]
vertices = Base.unsafe_wrap(Matrix{Float64},pointer(cp.vertices),(3,Int(cp.num_vertices[])))
fp.local_num_quadrants[]