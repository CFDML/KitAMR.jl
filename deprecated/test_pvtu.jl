using WriteVTK,MPI

MPI.Init()
if MPI.Comm_rank(MPI.COMM_WORLD) == 0
    vertices = [0.0 0.0;2.0 0.0;0.0 2.0;2.0 2.0;0.0 1.0;2.0 1.0]'
    ids = Vector{Vector{Int}}(undef,3)
    ids[1] = [1,12,14,10]
    ids[2] = [12,7,11,14]
    ids[3] = [10,14,13,5]
    cells = [MeshCell(PolyData.Polys(),ids[i]) for i in eachindex(ids)]
    pvtk_grid("test_polygon",vertices,cells;part = MPI.Comm_rank(MPI.COMM_WORLD)+1,
    nparts = MPI.Comm_size(MPI.COMM_WORLD),ismain = (MPI.Comm_rank(MPI.COMM_WORLD)==0),ghost_level = 1) do vtk
        vtk["test",VTKCellData()] = collect(1:3)
    end
end
if MPI.Comm_rank(MPI.COMM_WORLD) == 1
    vertices = [1.0 0.0;1.0 2.0;1.0 1.0;0.0 0.5;1.0 0.5;0.5 0.0;0.5 1.0;0.5 0.5]'
    ids = Vector{Vector{Int}}(undef,4)
    ids[1] = [14,11,9,13]
    ids[2] = [7,2,6,9,11]
    ids[3] = [5,13,9,8,3]
    ids[4] = [9,6,4,8]
    cells = [MeshCell(PolyData.Polys(),ids[i]) for i in eachindex(ids)]
    pvtk_grid("test_polygon",vertices,cells;part = MPI.Comm_rank(MPI.COMM_WORLD)+1,
    nparts = MPI.Comm_size(MPI.COMM_WORLD),ismain = (MPI.Comm_rank(MPI.COMM_WORLD)==0),ghost_level = 1) do vtk
        vtk["test",VTKCellData()] = collect(1:4)
    end
end
# vertices = [0.0 0.0;2.0 0.0;0.0 2.0;2.0 2.0;0.0 1.0;2.0 1.0;1.0 0.0;1.0 2.0;1.0 1.0;0.0 0.5;1.0 0.5;0.5 0.0;0.5 1.0;0.5 0.5]'
# ids = Vector{Vector{Int}}(undef,7)
# ids[1] = [1,12,14,10]
# ids[2] = [12,7,11,14]
# ids[3] = [10,14,13,5]
# ids[4] = [14,11,9,13]
# ids[5] = [7,2,6,9,11]
# ids[6] = [5,13,9,8,3]
# ids[7] = [9,6,4,8]
# cells = [MeshCell(PolyData.Polys(),ids[i]) for i in 1:7]
# vtk_grid("test_polygon",vertices,cells;append=false) do vtk
#     vtk["test",VTKCellData()] = collect(1:7)
# end