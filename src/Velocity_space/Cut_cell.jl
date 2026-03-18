function gaussian_area(A::AbstractMatrix) # DIMxN
    area = 0
    for i in axes(A,2)
        j = i%size(A,2)+1
        area+=@views det(A[:,[i,j]])
    end
    return 0.5*abs(area)
end

"""
2D cut cell.
"""
function cut_rect(n::Vector{Float64},vertices::Vector{Vector{Float64}}) # The splitting line at the boundary with normal direction n. Return the areas of the two part. The former one locates at the direction n.
    # The cases where n is oriented to along the axes are excluded in advance.
    #=
    y
#     7(8)|-----6-----|5(6)
#         |     |     |
#         8-----|-----4
#         |     |     |
#     1(2)|-----2-----|3(4) x
    =#
    points = Vector{Float64}[];indices = Int[];i = 0
    xmin,xmax,ymin,ymax = vertices[1][1],vertices[2][1],vertices[1][2],vertices[3][2]
    y = -n[1]/n[2]*xmin
    if y<ymax&&y>ymin
        push!(points,[xmin,y])
        push!(indices,8)
        i+=1
    end
    y = -n[1]/n[2]*xmax
    if y<ymax&&y>ymin
        push!(points,[xmax,y])
        push!(indices,4)
        i+=1
    end
    x = -n[2]/n[1]*ymin
    if x<xmax&&x>xmin
        push!(points,[x,ymin])
        push!(indices,2)
        i+=1
    end
    x = -n[2]/n[1]*ymax
    if x<xmax&&x>xmin
        push!(points,[x,ymax])
        push!(indices,6)
        i+=1
    end
    clp = findall(x->dot(x,n)==0.,vertices) # Number of vertices lying on the discontinuity
    if length(clp)+i<2
        return false,0.,0.
    else
        append!(points,vertices)
        append!(indices,CLP) # CLP is defined in dim.jl
        sid = sortperm(indices)
        indices = indices[sid]
        points = points[sid]
        aid = findfirst(x->iseven(x)||in(x,CLP[clp]),indices)
        bid = findlast(x->iseven(x)||in(x,CLP[clp]),indices)
        A = Matrix{Float64}(undef,2,bid-aid+1);B = Matrix{Float64}(undef,2,length(indices)-(bid-aid)+1)
        for i in aid:bid
            A[:,i-aid+1] .= points[i]
        end
        for i in 0:length(indices)-(bid-aid)
            index = bid+i>length(indices) ? (bid+i)%length(indices) : bid+i
            B[:,i+1] .= points[index]
        end
        l = points[bid]-points[aid]
        if n[1]*l[2]-n[2]*l[1]<0
            gas_weight = gaussian_area(A)
            return true,gas_weight,((xmax-xmin)*(ymax-ymin))-gas_weight # solid first
        else
            solid_weight = gaussian_area(A)
            return true,((xmax-xmin)*(ymax-ymin))-solid_weight,solid_weight
        end
    end
end

"""
3D cut cell.
"""
function cut_cube_rotate(n::Vector{Float64})
    C = zeros(2,3)
    e1 = [1.,0.,0.]
    e1 = cross(n,e1);e1./=norm(e1)
    e2 = cross(n,e1);e2./=norm(e2)
    if dot(cross(e1,e2),n)>0
        C[1,:].=e1;C[2,:].=e2
    else
        C[2,:].=e1;C[1,:].=e2
    end
    return C
end

"""
Calculate the flux of the vector (x,y,z)/3 through the vertical faces composed by intersecting points and vertices.
Normally, the cases where the connecting line of the two neighboring intersecting pionts are parallel to grid lines should have been got rid of before the function is called.
"""
function Vertical_Volume_Flux(points::AbstractVector{Vector{Float64}},midpoint::Vector{Float64},vertices::AbstractMatrix{Float64})
    H = 0.;N = length(points)
    for i in eachindex(points)
        p1 = points[i%N+1];p2 = points[i]
        dx= p1-p2
        id = findall(x->abs(x)<3.0*eps(),dx)
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
function vertices_sweep!(midpoint,ddu,vertices) # Clean the eps in vertices.
    if any(i->abs(midpoint[i])≈0.5*ddu[i],1:length(midpoint))
        for i in eachindex(vertices)
            abs(vertices[i]) < EPS && (vertices[i] = 0.)
        end
    end
    return nothing
end
function cut_cube(n::Vector{Float64},C::Matrix{Float64},midpoint::Vector{Float64},ddu::Vector{Float64},vertices::Matrix{Float64}) # 3.627 μs (215 allocations: 8.45 KiB). Acceptable?
    vertices_sweep!(midpoint,ddu,vertices)
    vltable = [[1,2],[3,4],[7,8],[5,6],[1,3],[2,4],[6,8],[5,7],[1,5],[2,6],[4,8],[3,7]] # vertices-edges table
    points = Vector{Vector{Float64}}(undef,6);index = 1
    dirs = permutedims(vertices)*n # what if vertices[i]=0.?
    for i in eachindex(vltable)
        d1 = dirs[vltable[i][1]]; d2 = dirs[vltable[i][2]]
        flag = d1*d2 # flag==0: cut any end of the edge; flag<0: cut the edge; flag>0: not cut the edge
        if min(abs(d1),abs(d2))<3.0*eps()
            if cld(i,4)==1 # avoid redundancy
                if abs(dirs[vltable[i][1]])<3.0*eps() # end A intersects
                    points[index] = vertices[:,vltable[i][1]];index+=1
                else # end B intersects
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
    index<4&&return false,0.,0.
    points = points[1:index-1]
    posid = findall(x->x>3.0*eps(),dirs)
    pos = length(posid)<4 ? true : false # H represents solid?
    if any(x->abs(x)<EPS,n) # Simple case
        dir = findfirst(x->abs(x)<EPS,n) 
        pid = findall(x->abs(x[dir]-vertices[dir])<EPS,points) # vertices[dir]: the dir-th component of the first vertex
        if pos
            vid = findall(x->abs(x[dir]-vertices[dir])<EPS&&dot(x,n)>EPS,eachcol(vertices)) # all vertices share the same face with the first one, but not cut by the boundary face
        else
            vid = findall(x->abs(x[dir]-vertices[dir])<EPS&&dot(x,n)<-EPS,eachcol(vertices)) # all vertices share the same face with the first one, but not cut by the boundary face
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
        if length(posid)==4 # isolated vertical face
            centers = @views sum(vertices[:,posid],dims = 2)./4.0
            dir = findfirst(i->abs(centers[i]-midpoint[i])≈0.5*ddu[i],1:3)
            if isnothing(dir)
                H = 0.
            else
                id1 = dir%3+1;id2 = (dir+1)%3+1
                H = pos ? sign(centers[dir]-midpoint[dir])*centers[dir]*ddu[id1]*ddu[id2]/3.0 : -sign(centers[dir]-midpoint[dir])*(2.0*midpoint[dir]-centers[dir])*ddu[id1]*ddu[id2]/3.0
            end
        else
            H = 0.
        end
        if pos
            @views H+=Vertical_Volume_Flux(points[id],midpoint,vertices[:,posid])
            return true,8*prod(midpoint-@view(vertices[:,1]))-H,H # gas first
        else
            negid = findall(x->x<-3.0*eps(),dirs)
            @views H+=Vertical_Volume_Flux(points[id],midpoint,vertices[:,negid])
            return true,H,8*prod(midpoint-@view(vertices[:,1]))-H
        end
    end
end