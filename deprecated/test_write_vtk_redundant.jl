using KitAMR,MPI,CairoMakie

MPI.Init()
config = KitAMR.read_config("configure.txt")
ps4est,amr = KitAMR.init(config);
global_data = amr.global_data



xmin,xmax,ymin,ymax,zmin,zmax = global_data.config.quadrature
Nx,Ny,Nz = global_data.config.vs_trees_num
AMR_VS_MAXLEVEL = global_data.config.solver.AMR_VS_MAXLEVEL

dx = (xmax - xmin) / Nx/2^AMR_VS_MAXLEVEL
dy = (ymax - ymin) / Ny/2^AMR_VS_MAXLEVEL
dz = (zmax - zmin) / Nz/2^AMR_VS_MAXLEVEL
D = [dx,dy,dz]

vs_data = amr.field.trees.data[end][end].vs_data
vertices = Matrix{Float64}(undef,3,8*length(vs_data.level))
midpoint = vs_data.midpoint
level = vs_data.level
dlevel = -(vs_data.level .- AMR_VS_MAXLEVEL)
cells = Vector{MeshCell}(undef,length(level))
for i in eachindex(level)
    for j in 1:8
        @. vertices[:,(i-1)*8+j] = midpoint[i,:]+RMT[3][j]*2^dlevel[i]*D/2
    end
    cells[i] = MeshCell(VTKCellTypes.VTK_VOXEL,(1:8).+8*(i-1))
end


vtk_grid("test_write_vtk",vertices,cells;append=false) do vtk
    vtk["DistrubutionFunction",VTKCellData()] = vs_data.df[:,1]
    vtk["level",VTKCellData()] = level
end