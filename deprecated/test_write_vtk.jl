using KitAMR,MPI,CairoMakie,WriteVTK

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

vertices_C = Array{Vector{Float64}}(undef, 2^AMR_VS_MAXLEVEL*Nx + 1, 2^AMR_VS_MAXLEVEL*Ny + 1, 2^AMR_VS_MAXLEVEL*Nz + 1)
for k = 1:2^AMR_VS_MAXLEVEL*Nz+1
    for j = 1:2^AMR_VS_MAXLEVEL*Ny+1
        for i = 1:2^AMR_VS_MAXLEVEL*Nx+1
            vertices_C[i, j, k] = [xmin + (i - 1) * dx, ymin + (j - 1) * dy, zmin + (k - 1) * dz]
        end
    end
end

vL = LinearIndices(vertices_C)

vertices = reshape(vertices_C, :)
V = Matrix{Float64}(undef,3,length(vertices))
for i in eachindex(vertices)
    V[:,i].=vertices[i]
end

vs_data = amr.field.trees.data[end][end].vs_data
pm = [[-1,-1,-1],[1,-1,-1],[-1,1,-1],[1,1,-1],[-1,-1,1],[1,-1,1],[-1,1,1],[1,1,1]]
midpoint = vs_data.midpoint
level = vs_data.level
dlevel = -(vs_data.level .- AMR_VS_MAXLEVEL)
v = Vector{Int}(undef,8)
# function make_polyhedron(v)
#      VTKPolyhedron(v,(v[1],v[4],v[3],v[2]),(v[1],v[5],v[8],v[4]),(v[5],v[6],v[7],v[8]),(v[6],v[2],v[3],v[7]),(v[1],v[2],v[6],v[5]),(v[3],v[4],v[8],v[7]))
# end
# cells = Vector{VTKPolyhedron}(undef,length(level))
# for i in eachindex(level)
#     for j in eachindex(v)
#         vinner = Vector{Int}(undef,3)
#         for k in 1:3
#             vinner[k] = round((midpoint[i,k]+pm[j][k]*2^dlevel[i]*D[k]/2-global_data.config.quadrature[2*k-1])/D[k])+1
#         end
#         v[j] = vL[vinner...]
#     end
#     cells[i] = make_polyhedron(copy(v))
# end

cells = Vector{MeshCell}(undef,length(level))
vinner = Vector{Int}(undef,3)
for i in eachindex(level)
    for j in eachindex(v)
        for k in 1:3
            vinner[k] = round((midpoint[i,k]+pm[j][k]*2^dlevel[i]*D[k]/2-global_data.config.quadrature[2*k-1])/D[k])+1
        end
        v[j] = vL[vinner...]
    end
    cells[i] = MeshCell(VTKCellTypes.VTK_VOXEL,copy(v))
end

vtk_grid("test_write_vtk",V,cells;append=false) do vtk
    vtk["DistrubutionFunction",VTKCellData()] = vs_data.df[:,1]
    vtk["level",VTKCellData()] = level
end

points = permutedims(Matrix{Float64}([0 0 0;1 0 0;0 1 0;1 1 0;0 0 1;1 0 1;0 1 1;1 1 1]))
cells = MeshCell(VTKCellTypes.VTK_VOXEL,1:8)
vtk_grid("test_write_vtk",points,[cells];append=false) do vtk
end