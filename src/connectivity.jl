function Cartesian_connectivity(Nx, Ny, xmin, xmax, ymin, ymax)
    vertices_C = Array{NTuple{2,Float64}}(undef, Nx + 1, Ny + 1)
    dx = (xmax - xmin) / Nx
    dy = (ymax - ymin) / Ny
    for j = 1:(Ny+1)
        for i = 1:(Nx+1)
            vertices_C[i, j] = (xmin + (i - 1) * dx, ymin + (j - 1) * dy)
        end
    end
    vL = LinearIndices(vertices_C)
    vertices = reshape(vertices_C, :)
    cell2ver = Array{NTuple{4,Int32}}(undef, Nx, Ny)
    for j in axes(cell2ver, 2)
        for i in axes(cell2ver, 1)
            cell2ver[i, j] = (vL[i, j], vL[i+1, j], vL[i, j+1], vL[i+1, j+1])
        end
    end
    cells = reshape(cell2ver, :)
    connectivity = GC.@preserve vertices cells Connectivity{4}(vertices, cells)
end
