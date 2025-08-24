function Cartesian_connectivity(Nx, Ny, xmin, xmax, ymin, ymax)
    vertices_C = Array{NTuple{2,Float64}}(undef, Nx + 1, Ny + 1)
    dx = (xmax - xmin) / Nx
    dy = (ymax - ymin) / Ny
    for j = 1:Ny+1
        for i = 1:Nx+1
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
function Cartesian_connectivity(Nx,Ny,Nz,xmin,xmax,ymin,ymax,zmin,zmax)
    vertices_C = Array{NTuple{3,Float64}}(undef, Nx + 1, Ny + 1, Nz + 1)
    dx = (xmax - xmin) / Nx
    dy = (ymax - ymin) / Ny
    dz = (zmax - zmin) / Nz
    for k = 1:Nz+1
        for j = 1:Ny+1
            for i = 1:Nx+1
                vertices_C[i, j, k] = (xmin + (i - 1) * dx, ymin + (j - 1) * dy, zmin + (k - 1) * dz)
            end
        end
    end
    vL = LinearIndices(vertices_C)
    vertices = reshape(vertices_C, :)
    cell2ver = Array{NTuple{8,Int32}}(undef, Nx, Ny, Nz)
    for k in axes(cell2ver, 3)
        for j in axes(cell2ver, 2)
            for i in axes(cell2ver, 1)
                cell2ver[i, j, k] = (vL[i, j, k], vL[i+1, j, k], vL[i, j+1, k], vL[i+1, j+1, k],
                                     vL[i, j, k+1], vL[i+1, j, k+1], vL[i, j+1, k+1], vL[i+1, j+1, k+1])
                # cell2ver[i, j, k] = (vL[i, j, k], vL[i, j, k+1], vL[i+1, j, k], vL[i+1, j, k+1],
                #                      vL[i, j+1, k+1], vL[i, j+1, k+1], vL[i+1, j+1, k], vL[i+1, j+1, k+1])
            end
        end
    end
    cells = reshape(cell2ver, :)
    connectivity = GC.@preserve vertices cells Connectivity{8}(vertices, cells)
end
function set_connectivity(global_data::Global_Data{DIM}) where{DIM}
    domain = global_data.config.domain

    periodic_dirs = ntuple(i -> begin
        first_boundary = (i-1)*2 + 1  # 1, 3, 5 for i = 1, 2, 3
        second_boundary = first_boundary + 1  # 3, 5, 7 for i = 1, 2, 3
        nameof(typeof(domain[first_boundary]).parameters[1]) == :Period &&
        nameof(typeof(domain[second_boundary]).parameters[1]) == :Period
    end, DIM)
    
    if periodic_dirs != ntuple(_ -> false, DIM)
        connectivity_ps = set_periodic_connectivity(global_data, periodic_dirs)
    else
        connectivity_ps = Cartesian_connectivity(global_data.config.trees_num..., global_data.config.geometry...)
    end
    return connectivity_ps
end

import ..P4estTypes.unsafe_vertices
function set_periodic_connectivity(global_data::Global_Data{DIM}, periodic_dirs) where{DIM}
    geometry = global_data.config.geometry
    bounds = if DIM == 2
        ntuple(i -> begin
            if i <= 2
                min_idx = (i-1)*2 + 1  # 1, 3 for i = 1, 2
                max_idx = min_idx + 1   # 2, 4 for i = 1, 2
                (geometry[min_idx], geometry[max_idx])  # (xmin,xmax), (ymin,ymax)
            else
                (0.0, 0.0)
            end
        end, 3)
    else
        ntuple(i -> begin
            min_idx = (i-1)*2 + 1  # 1, 3, 5 for i = 1, 2, 3
            max_idx = min_idx + 1   # 2, 4, 6 for i = 1, 2, 3
            (geometry[min_idx], geometry[max_idx])  # (xmin,xmax), (ymin,ymax), (zmin,zmax)
        end, 3)
    end
    
    connectivity_ps = P4estTypes.brick(
        Tuple(global_data.config.trees_num),
        periodic_dirs
    )
    vertices = unsafe_vertices(connectivity_ps)
    for i in eachindex(vertices)
        normalized_coords = vertices[i]
        
        vertices[i] = ntuple(j -> begin
            if j <= DIM
                normalized_coords[j] * (bounds[j][2] - bounds[j][1])/global_data.config.trees_num[j] + bounds[j][1]
            else
                0.0
            end
        end, 3)
    end
    return connectivity_ps
end