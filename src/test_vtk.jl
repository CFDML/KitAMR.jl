using KitAMR,MPI,CairoMakie

MPI.Init()
config = KitAMR.read_config("configure.txt")
ps4est,amr = KitAMR.init(config);
p4est_geo = p8est_geometry_new_connectivity(pointer(PointerWrapper(ps4est).connectivity))
cont = p8est_vtk_context_new(ps4est,"testZ")
p8est_vtk_context_set_geom(cont,p4est_geo)
p8est_vtk_context_set_continuous(cont,1)
cont = p8est_vtk_write_header(cont)
fp = PointerWrapper(ps4est)
num_quads = fp.local_num_quadrants[]
p = sc_malloc(-1,num_quads*8)
pa = Base.unsafe_wrap(Vector{Cdouble},Ptr{Cdouble}(p),num_quads)
ppa = Base.pointer_from_objref(pa)
psc = sc_array_new_data(p,8,num_quads)
function init_cell_data_kernel(ip,data,dp)
    pa = Base.unsafe_pointer_to_objref(data)
    qid = KitAMR.global_quadid(ip)+1
    ps_data = Base.unsafe_pointer_to_objref(pointer(dp.ps_data))
    pa[qid] = ps_data.midpoint[3]
end
function init_cell_data(info,data)
    KitAMR.AMR_volume_iterate(info, data, KitAMR.P4est_PS_Data, init_cell_data_kernel)
end
GC.@preserve pa KitAMR.AMR_4est_volume_iterate(ps4est,amr.global_data.forest.ghost,ppa,init_cell_data)
cont = p8est_vtk_write_cell_dataf(cont,1,1,1,0,1,0,"testZ",psc,cont)
retval = p8est_vtk_write_footer(cont)
