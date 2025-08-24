using LinearAlgebra, Statistics
function gaussian_area(A::AbstractMatrix) # DIMxN
    area = 0
    for i in axes(A,2)
        j = i%size(A,2)+1
        area+=@views det(A[:,[i,j]])
    end
    return 0.5*abs(area)
end
function cut_cube_CA(n::Vector{Float64})
    A = zeros(2);C = zeros(2,3)
    e1 = [1.,0.,0.]
    e1 = cross(n,e1);e1./=norm(e1)
    e2 = cross(n,e1);e2./=norm(e2)
    if dot(cross(e1,e2),n)>0
        C[1,:].=e1;C[2,:].=e2
        A[1]=dot(e1,n);A[2]=dot(e2,n)
    else
        C[2,:].=e1;C[1,:].=e2
        A[2]=dot(e1,n);A[1]=dot(e2,n)
    end
    return C,A/3.0
end
function Volume_Flux(points::AbstractMatrix{Float64},A::Vector{Float64}) # Calculate the flux (x,y,z)/3 through arbitrary polygon to evaluate the volume of an arbitrary polyhedron.
    # points: The sorted vertices of the polygon.
    # A: The preconfigured coefficient
    H = 0.
    N = size(points,2)
    for i in axes(points,2)
        dx,dy = @views points[:,(i)%N+1]-points[:,i]
        x,y = @views points[:,i]
        C0 = x*dy*(0.5*A[1]*x+A[2]*y)
        C1 = dy*(x*dy*A[2]+dx*(A[1]*x+A[2]*y))
        C2 = dx*dy*(0.5*dx*A[1]+dy*A[2])
        H+=C0+0.5*C1+C2/3.0
    end
    return H
end
function Vertical_Volume_Flux(points::AbstractVector{Vector{Float64}},midpoint::Vector{Float64},vertices::AbstractMatrix{Float64})
    H = 0.;N = length(points)
    for i in eachindex(points)
        p1 = points[i%N+1];p2 = points[i]
        dx= p1-p2
        EPS = 1e-12
        id = findall(x->abs(x)<EPS,dx)
        if length(id)>1
            throw(`Cut-cube Error!`)
        else
            dir = id[1]
            vid = findall(x->abs(x[dir]-points[i][dir])<EPS,eachcol(vertices))
            A = hcat(vertices[:,vid],p1,p2)
            center = vec(mean(A,dims=2))
            id1 = dir%3+1;id2 = (dir+1)%3+1
            phi=@views [atan(x[id2]-center[id2],x[id1]-center[id1]) for x in eachcol(A)]
            id_gauss = sortperm(phi)
            H+=sign(p1[dir]-midpoint[dir])*p1[dir]/3.0*gaussian_area(@views A[[id1,id2],id_gauss])
        end
    end
    return H
end
function cut_cube(n::Vector{Float64},C::Matrix{Float64},A::Vector{Float64},midpoint::Vector{Float64},vertices::Matrix{Float64}) # 3.627 μs (215 allocations: 8.45 KiB). Acceptable?
    vltable = [[1,2],[3,4],[7,8],[5,6],[1,3],[2,4],[6,8],[5,7],[1,5],[2,6],[4,8],[3,7]] # vertices-edges table
    points = Vector{Vector{Float64}}(undef,6);index = 1
    dirs = permutedims(vertices)*n
    EPS = 1e-12
    for i in eachindex(vltable)
        flag = dirs[vltable[i][1]]*dirs[vltable[i][2]]
        if abs(flag)<EPS
            if cld(i,4)==1 # avoid redundancy
                if abs(dirs[vltable[i][1]])<EPS # end A intersects
                    points[index] = vertices[:,vltable[i][1]];index+=1
                end
                if abs(dirs[vltable[i][2]])<EPS # end B intersects
                    points[index] = vertices[:,vltable[i][2]];index+=1
                end
            end
        elseif flag<0 # intersects between the two ends
            point = vertices[:,vltable[i][1]]
            dir = cld(i,4);point[dir]=0.
            point[dir] = -dot(point,n)/(n[dir])
            points[index] = point;index+=1
        end
    end
    @show points
    index<4&&return false,0.,0.
    points = points[1:index-1]
    posid = findall(x->x>EPS,dirs)
    pos = length(posid)<4 ? true : false
    if any(x->abs(x)<EPS,n) # Simple case
        dir = findfirst(x->abs(x)<EPS,n) 
        pid = findall(x->abs(x[dir]-vertices[dir])<EPS,points)
        if pos
            vid = findall(x->abs(x[dir]-vertices[dir])<EPS&&dot(x,n)>EPS,eachcol(vertices))
        else
            vid = findall(x->abs(x[dir]-vertices[dir])<EPS&&dot(x,n)<-EPS,eachcol(vertices))
        end
        A = hcat(vertices[:,vid],points[pid]...)
        center = vec(mean(A,dims=2))
        id1 = dir%3+1;id2 = (dir+1)%3+1
        phi=@views [atan(x[id2]-center[id2],x[id1]-center[id1]) for x in eachcol(A)]
        id_gauss = sortperm(phi)
        H = 2*(midpoint[dir]-vertices[dir,1])*gaussian_area(@views A[[id1,id2],id_gauss])
        if pos
            return true,8*prod(midpoint-@view(vertices[:,1]))-H,H # gas first
        else
            return true,H,8*prod(midpoint-@view(vertices[:,1]))-H
        end
    else
        local_points = Matrix{Float64}(undef,2,length(points))
        for i in eachindex(points)
            local_points[:,i].=C*points[i]
        end
        center = vec(mean(local_points,dims=2))
        phi=@views [atan(x[2]-center[2],x[1]-center[1]) for x in eachcol(local_points)]
        id = sortperm(phi)
        @views H = pos ? Volume_Flux(local_points[:,id],A) : -Volume_Flux(local_points[:,id],A)
        if pos
            # @show vertices[:,posid]
            @views H+=Vertical_Volume_Flux(points[id],midpoint,vertices[:,posid])
            # if posid==4 # The case that the vertical face is purely composed by vertices
            #     @views v1 = vertices[:,posid[1]]
            #     @views v2 = vertices[:,posid[2]]
            #     @views v3 = vertices[:,posid[3]]
            #     pface_dir = findfirst(x->abs(x)>EPS,cross(v1-v2,v1-v3))
            #     H+=vertices[pface_dir,posid[1]]/3.0*gaussian_area(@views vertices[:,posid])
            # end
            return true,8*prod(midpoint-@view(vertices[:,1]))-H,H # gas first
        else
            negid = findall(x->x<-EPS,dirs)
            @views H+=Vertical_Volume_Flux(points[id],midpoint,vertices[:,negid])
            return true,H,8*prod(midpoint-@view(vertices[:,1]))-H
        end
    end
end
midpoint = [0.5,0.5,0.5];n = [-2/sqrt(5),1/sqrt(5),0]
vertices = [0. 0 0;1 0 0;0 1 0;1 1 0;0 0 1;1 0 1;0 1 1;1 1 1] |> permutedims
C,A = cut_cube_CA(n)
flag,gas_weight,solid_weight = cut_cube(n,C,A,midpoint,vertices)

midpoint = [0.5,1.5,0.5];n = [-2,0.5,0]; n ./= norm(n)
vertices = [0. 0 0;1 0 0;0 1 0;1 1 0;0 0 1;1 0 1;0 1 1;1 1 1]
vertices[:,2].+=1.0
vertices = permutedims(vertices)
C,A = cut_cube_CA(n)
flag,gas_weight,solid_weight = cut_cube(n,C,A,midpoint,vertices)

midpoint = [0.5,0.5,0.5];n = [-1,1,1.]; n ./= norm(n)
vertices = [0. 0 0;1 0 0;0 1 0;1 1 0;0 0 1;1 0 1;0 1 1;1 1 1] |> permutedims
C,A = cut_cube_CA(n)
flag,gas_weight,solid_weight = cut_cube(n,C,A,midpoint,vertices)