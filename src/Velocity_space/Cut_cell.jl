function gaussian_area(A::AbstractMatrix) # DIMxN
    area = 0
    for i in axes(A,2)
        j = i%size(A,2)+1
        area+=@views det(A[:,[i,j]])
    end
    return 0.5*abs(area)
end

"""
2D cut cell. `gas_weight` is the area in `dot(v,n)<0`; `solid_weight` is
the area in `dot(v,n)>0`.
"""
function _same_cut_point(x1::Float64,y1::Float64,x2::Float64,y2::Float64,tol::Float64)
    return abs(x1-x2) <= tol && abs(y1-y2) <= tol
end
function _append_cut_point(count::Int,area2::Float64,first_x::Float64,first_y::Float64,last_x::Float64,last_y::Float64,x::Float64,y::Float64,tol::Float64)
    if count == 0
        count = 1
        first_x = x;first_y = y
        last_x = x;last_y = y
    elseif !_same_cut_point(last_x,last_y,x,y,tol)
        area2 += last_x*y-last_y*x
        count += 1
        last_x = x;last_y = y
    end
    return count,area2,first_x,first_y,last_x,last_y
end
function _cut_rect_tol(n::Vector{Float64},vertices::Vector{Vector{Float64}})
    scale = 1.0
    for vertex in vertices
        for x in vertex
            scale = max(scale,abs(x))
        end
    end
    return max(EPS,128.0*eps(Float64)*max(sqrt(n[1]*n[1]+n[2]*n[2]),1.0)*scale)
end
function _rect_area(vertices::Vector{Vector{Float64}})
    area = 0.0
    @inbounds for i in 1:4
        j = i == 4 ? 1 : i+1
        area += vertices[i][1]*vertices[j][2]-vertices[i][2]*vertices[j][1]
    end
    return 0.5*abs(area)
end
function _clip_rect_negative_halfplane_area(n::Vector{Float64},vertices::Vector{Vector{Float64}},tol::Float64)
    count = 0
    area2 = 0.0
    first_x = 0.0;first_y = 0.0
    last_x = 0.0;last_y = 0.0
    n1 = n[1];n2 = n[2]
    prev = vertices[4]
    prev_x = prev[1];prev_y = prev[2]
    dprev = n1*prev_x+n2*prev_y
    prev_inside = dprev <= tol
    @inbounds for i in 1:4
        curr = vertices[i]
        curr_x = curr[1];curr_y = curr[2]
        dcurr = n1*curr_x+n2*curr_y
        curr_inside = dcurr <= tol
        if curr_inside != prev_inside
            denom = dprev-dcurr
            if abs(denom) > tol
                theta = clamp(dprev/denom,0.0,1.0)
                count,area2,first_x,first_y,last_x,last_y = _append_cut_point(
                    count,area2,first_x,first_y,last_x,last_y,
                    prev_x+theta*(curr_x-prev_x),prev_y+theta*(curr_y-prev_y),tol
                )
            end
        end
        if curr_inside
            count,area2,first_x,first_y,last_x,last_y = _append_cut_point(
                count,area2,first_x,first_y,last_x,last_y,curr_x,curr_y,tol
            )
        end
        prev_x = curr_x;prev_y = curr_y
        dprev = dcurr
        prev_inside = curr_inside
    end
    if count>1&&!_same_cut_point(first_x,first_y,last_x,last_y,tol)
        area2 += last_x*first_y-last_y*first_x
    end
    return count,0.5*abs(area2)
end
function cut_rect(n::Vector{Float64},vertices::Vector{Vector{Float64}})
    total_weight = _rect_area(vertices)
    total_weight <= 0.0&&return false,0.,0.
    tol = _cut_rect_tol(n,vertices)
    count,gas_weight = _clip_rect_negative_halfplane_area(n,vertices,tol)
    count<3&&return false,0.,0.
    gas_weight = clamp(gas_weight,0.0,total_weight)
    solid_weight = total_weight-gas_weight
    area_tol = max(tol*tol,128.0*eps(Float64)*max(total_weight,1.0))
    if gas_weight <= area_tol||solid_weight <= area_tol
        return false,0.,0.
    end
    return true,gas_weight,solid_weight
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
        dx= abs.(p1-p2)
        _,dir = findmin(dx)
        if dx[dir]>1e-6
            throw(`Cut-cube Error!`)
        end
        vid = findall(x->abs(x[dir]-points[i][dir])<EPS,eachcol(vertices))
        A = hcat(vertices[:,vid],p1,p2)
        center = vec(mean(A,dims=2))
        id1 = dir%3+1;id2 = (dir+1)%3+1
        phi=@views [atan(x[id2]-center[id2],x[id1]-center[id1]) for x in eachcol(A)]
        id_gauss = sortperm(phi)
        H+=sign(p1[dir]-midpoint[dir])*p1[dir]/3.0*gaussian_area(@views A[[id1,id2],id_gauss])
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
const CUT_CUBE_EDGE_VERTICES = ((1,2),(3,4),(7,8),(5,6),(1,3),(2,4),(6,8),(5,7),(1,5),(2,6),(4,8),(3,7))
function cut_cube(n::Vector{Float64},C::Matrix{Float64},midpoint::Vector{Float64},ddu::Vector{Float64},vertices::Matrix{Float64}) # 3.627 μs (215 allocations: 8.45 KiB). Acceptable?
    vertices_sweep!(midpoint,ddu,vertices)
    points = Vector{Vector{Float64}}(undef,6);index = 1
    dirs = MVector{8,Float64}(undef)
    @inbounds for i in 1:8
        dirs[i] = vertices[1,i]*n[1]+vertices[2,i]*n[2]+vertices[3,i]*n[3]
    end
    for i in eachindex(CUT_CUBE_EDGE_VERTICES)
        edge = CUT_CUBE_EDGE_VERTICES[i]
        d1 = dirs[edge[1]]; d2 = dirs[edge[2]]
        flag = d1*d2 # flag==0: cut any end of the edge; flag<0: cut the edge; flag>0: not cut the edge
        if min(abs(d1),abs(d2))<3.0*eps()
            if cld(i,4)==1 # avoid redundancy
                if abs(dirs[edge[1]])<3.0*eps() # end A intersects
                    points[index] = vertices[:,edge[1]];index+=1
                else # end B intersects
                    points[index] = vertices[:,edge[2]];index+=1
                end
            end
        elseif flag<0 # intersects between the two ends
            point = vertices[:,edge[1]]
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
