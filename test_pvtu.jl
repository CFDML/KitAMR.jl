using WriteVTK,MPI
MPI.Init()
if MPI.Comm_rank(MPI.COMM_WORLD)==0
    points = [[0.0,0.0],[0.5,0.0],[1.0,0.0],[0.0,0.5],[0.5,0.5],[1.0,0.5],[0.,-0.5],[0.5,-0.5],[1.0,-0.5]]
    vertices = zeros(2,length(points))
    for i in eachindex(points)
        vertices[:,i].=points[i]
    end
    vtk_cnn = [[1,2,5,4],[2,3,6,5],[7,8,2,1],[8,9,3,2]]
    cells = [MeshCell(PolyData.Polys(),cnn) for cnn in vtk_cnn]
    pvtk_grid("test0",vertices,cells;part = 1, nparts = 2, ghost_level = 1) do pvtk
        pvtk["test"] = [1.,2.,-1.,0.]
    end
else
    points = [[0.0,0.0],[0.5,0.0],[1.0,0.0],[0.0,0.5],[0.5,0.5],[1.0,0.5],[0.0,1.0],[0.5,1.0],[1.0,1.0]]
    vertices = zeros(2,length(points))
    for i in eachindex(points)
        vertices[:,i].=points[i]
    end
    vtk_cnn = [[1,2,5,4],[2,3,6,5],[4,5,8,7],[5,6,9,8]]
    cells = [MeshCell(PolyData.Polys(),cnn) for cnn in vtk_cnn]
    pvtk_grid("test0",vertices,cells;part = 2, nparts = 2, ghost_level = 1) do pvtk
        pvtk["test"] = [1.,2.,3.,4.]
    end
end
MPI.Finalize()

function test(;c::Float64,kwargs...)
    kwargs[:a]+kwargs[:b]+c
end
test(;a = 1.,b = 2., c= 3.)