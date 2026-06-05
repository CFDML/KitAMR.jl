using LinearAlgebra, Statistics
function cut_cube_rotate(n::Vector{Float64})
    C = zeros(2, 3)
    e1 = [1.0, 0.0, 0.0]
    e1 = cross(n, e1);
    e1./=norm(e1)
    e2 = cross(n, e1);
    e2./=norm(e2)
    if dot(cross(e1, e2), n)>0
        C[1, :].=e1;
        C[2, :].=e2
    else
        C[2, :].=e1;
        C[1, :].=e2
    end
    return C
end
function Vertical_Volume_Flux(
    points::AbstractVector{Vector{Float64}},
    midpoint::Vector{Float64},
    vertices::AbstractMatrix{Float64},
)
    H = 0.0;
    N = length(points)
    for i in eachindex(points)
        p1 = points[i%N+1];
        p2 = points[i]
        dx = p1-p2
        id = findall(x->abs(x)<3.0*eps(), dx)
        if length(id)>1
            throw(`Cut-cube Error!`)
        else
            dir = id[1]
            vid = findall(x->abs(x[dir]-points[i][dir])<EPS, eachcol(vertices))
            A = hcat(vertices[:, vid], p1, p2)
            center = vec(mean(A, dims = 2))
            id1 = dir%3+1;
            id2 = (dir+1)%3+1
            phi=@views [atan(x[id2]-center[id2], x[id1]-center[id1]) for x in eachcol(A)]
            id_gauss = sortperm(phi)
            H+=sign(p1[dir]-midpoint[dir])*p1[dir]/3.0*gaussian_area(
                @views A[[id1, id2], id_gauss]
            )
        end
    end
    return H
end
function vertices_sweep!(midpoint, ddu, vertices) # Clean the eps in vertices.
    if any(i->abs(midpoint[i])≈0.5*ddu[i], 1:length(midpoint))
        for i in eachindex(vertices)
            abs(vertices[i]) < EPS && (vertices[i] = 0.0)
        end
    end
    return nothing
end
function cut_cube(
    n::Vector{Float64},
    C::Matrix{Float64},
    midpoint::Vector{Float64},
    ddu::Vector{Float64},
    vertices::Matrix{Float64},
) # 3.627 μs (215 allocations: 8.45 KiB). Acceptable?
    vertices_sweep!(midpoint, ddu, vertices)
    vltable = [
        [1, 2],
        [3, 4],
        [7, 8],
        [5, 6],
        [1, 3],
        [2, 4],
        [6, 8],
        [5, 7],
        [1, 5],
        [2, 6],
        [4, 8],
        [3, 7],
    ] # vertices-edges table
    points = Vector{Vector{Float64}}(undef, 6);
    index = 1
    dirs = permutedims(vertices)*n # what if vertices[i]=0.?
    @show dirs permutedims(vertices)
    for i in eachindex(vltable)
        d1 = dirs[vltable[i][1]];
        d2 = dirs[vltable[i][2]]
        flag = d1*d2 # flag==0: cut any end of the edge; flag<0: cut the edge; flag>0: not cut the edge
        if min(abs(d1), abs(d2))<3.0*eps()
            if cld(i, 4)==1 # avoid redundancy
                if abs(dirs[vltable[i][1]])<3.0*eps() # end A intersects
                    points[index] = vertices[:, vltable[i][1]];
                    index+=1
                else # end B intersects
                    points[index] = vertices[:, vltable[i][2]];
                    index+=1
                end
            end
        elseif flag<0 # intersects between the two ends
            point = vertices[:, vltable[i][1]]
            dir = cld(i, 4);
            point[dir]=0.0
            point[dir] = -dot(point, n)/(n[dir])
            points[index] = point;
            index+=1
        end
    end
    index<4&&return false, 0.0, 0.0
    points = points[1:(index-1)]
    posid = findall(x->x>3.0*eps(), dirs)
    pos = length(posid)<4 ? true : false # H represents solid?
    if any(x->abs(x)<EPS, n) # Simple case
        dir = findfirst(x->abs(x)<EPS, n)
        pid = findall(x->abs(x[dir]-vertices[dir])<EPS, points) # vertices[dir]: the dir-th component of the first vertex
        if pos
            vid =
                findall(x->abs(x[dir]-vertices[dir])<EPS&&dot(x, n)>EPS, eachcol(vertices)) # all vertices share the same face with the first one, but not cut by the boundary face
        else
            vid =
                findall(x->abs(x[dir]-vertices[dir])<EPS&&dot(x, n)<-EPS, eachcol(vertices)) # all vertices share the same face with the first one, but not cut by the boundary face
        end
        A = hcat(vertices[:, vid], points[pid]...)
        center = vec(mean(A, dims = 2))
        id1 = dir%3+1;
        id2 = (dir+1)%3+1
        phi=@views [atan(x[id2]-center[id2], x[id1]-center[id1]) for x in eachcol(A)]
        id_gauss = sortperm(phi)
        H = 2*(midpoint[dir]-vertices[dir, 1])*gaussian_area(@views A[[id1, id2], id_gauss])
        if pos
            return true, 8*prod(midpoint-@view(vertices[:, 1]))-H, H # gas first
        else
            return true, H, 8*prod(midpoint-@view(vertices[:, 1]))-H
        end
    else
        local_points = Matrix{Float64}(undef, 2, length(points))
        for i in eachindex(points)
            local_points[:, i].=C*points[i]
        end
        center = vec(mean(local_points, dims = 2))
        phi=@views [atan(x[2]-center[2], x[1]-center[1]) for x in eachcol(local_points)]
        id = sortperm(phi)
        if length(posid)==4 # isolated vertical face
            centers = @views sum(vertices[:, posid], dims = 2) ./ 4.0
            dir = findfirst(i->abs(centers[i]-midpoint[i])≈0.5*ddu[i], 1:3)
            if isnothing(dir)
                H = 0.0
            else
                id1 = dir%3+1;
                id2 = (dir+1)%3+1
                H =
                    pos ?
                    sign(centers[dir]-midpoint[dir])*centers[dir]*ddu[id1]*ddu[id2]/3.0 :
                    -sign(centers[dir]-midpoint[dir])*(2.0*midpoint[dir]-centers[dir])*ddu[id1]*ddu[id2]/3.0
            end
        else
            H = 0.0
        end
        if pos
            @show points[id]
            @views H+=Vertical_Volume_Flux(points[id], midpoint, vertices[:, posid])
            return true, 8*prod(midpoint-@view(vertices[:, 1]))-H, H # gas first
        else
            negid = findall(x->x<0.0, dirs)
            @views H+=Vertical_Volume_Flux(points[id], midpoint, vertices[:, negid])
            return true, H, 8*prod(midpoint-@view(vertices[:, 1]))-H
        end
    end
end
function gaussian_area(A::AbstractMatrix) # DIMxN
    area = 0
    for i in axes(A, 2)
        j = i%size(A, 2)+1
        area+=@views det(A[:, [i, j]])
    end
    return 0.5*abs(area)
end


function test_cut_cube(midpoint, ddu, n)
    ANTIVT = [
        [],
        [[-1, -1], [1, -1], [1, 1], [-1, 1]],
        [
            [-1, -1, -1],
            [1, -1, -1],
            [-1, 1, -1],
            [1, 1, -1],
            [-1, -1, 1],
            [1, -1, 1],
            [-1, 1, 1],
            [1, 1, 1],
        ],
    ] # anticlock vertices table
    C = cut_cube_rotate(n)
    vertices = zeros(3, 8)
    for j in axes(vertices, 2)
        vertices[:, j] = 0.5*ANTIVT[3][j] .* ddu+midpoint
    end
    cut_cube(n, C, midpoint, ddu, vertices)
end
#=
sphere bug
=#
# ddu = [0.91, 0.91, 0.91]
# midpoint = [-0.455, -6.825, -0.455]
# n = [-0.21875, 0.03125, 0.9752804083954522]

# n = [-0.294039855472166, -0.18628225226505013, -0.937464391795746]
# midpoint = [-0.22187499999999866, -0.221875, -0.221875]
# ddu = [0.44375, 0.44375, 0.44375]

n = [0.65625, -0.65625, 0.37238672774415577]
ddu = [0.91, 0.91, 0.91]
midpoint = [-6.825, -6.825, -0.455]
test_cut_cube(midpoint, ddu, n)
