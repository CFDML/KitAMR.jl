function pre_ps_refine_flag(boundary::AbstractCircle,midpoint::AbstractVector,ds::AbstractVector,kinfo::KInfo) # Circle type IB boundary flag
    boundary_flag(boundary,midpoint,ds,kinfo)
end
function boundary_flag(boundary::Circle,midpoint::AbstractVector,ds::AbstractVector,::KInfo) # Circle type IB boundary flag
    flag = 0
    for i = 1:4
        flag += norm(midpoint.+ds.*NMT[2][i].-boundary.center)>boundary.radius ? 1 : -1 # any neighbor cross boundary?
    end
    abs(flag)==4 && return false
    return true
end
function boundary_flag(boundary::Sphere,midpoint::AbstractVector,ds::AbstractVector,::KInfo) # Circle type IB boundary flag
    flag = 0
    for i = 1:6
        flag += norm(midpoint.+ds.*NMT[3][i].-boundary.center)>boundary.radius ? 1 : -1 # any neighbor cross boundary?
    end
    abs(flag)==6 && return false
    return true
end
function solid_flag(boundary::AbstractCircle,midpoint::AbstractVector) # Does midpoint locate at solid?
    return xor(norm(midpoint.-boundary.center)>boundary.radius,boundary.solid)
end
function solid_cell_flag(boundary::AbstractCircle,midpoint::AbstractVector,ds::AbstractVector,kinfo::KInfo,inside::Bool) # Ghost nodes, those are inside solid domain and immediately adjacent the boundary.
    (boundary_flag(boundary,midpoint,ds,kinfo) && inside) && return true
    return false
end
function calc_IB_ρw(aux_point::AbstractVector,bound::Circle,midpoint::AbstractMatrix,weight::AbstractVector,df::AbstractMatrix,vn::AbstractVector,Θ::AbstractVector)
    calc_IB_ρw_2D(aux_point,bound.bc,midpoint,weight,df,vn,Θ)
end
function calc_IB_ρw(aux_point::AbstractVector,bound::Sphere,midpoint::AbstractMatrix,weight::AbstractVector,df::AbstractMatrix,vn::AbstractVector,Θ::AbstractVector)
    calc_IB_ρw_3D(aux_point,bound.bc,midpoint,weight,df,vn,Θ)
end
function IB_prim(circle::AbstractCircle,aux_point::AbstractVector,ρw::Real)
    IB_prim(circle.bc,aux_point,ρw)
end
function calc_intersect(f_midpoint,s_midpoint,::Vector,::Int,circle::Circle)
    c = circle.center;r = circle.radius
    if abs(f_midpoint[1]-s_midpoint[1])<EPS
        t = acos((f_midpoint[1]-c[1])/r)
        if f_midpoint[2]>c[2]
            ap = [r*cos(t),r*sin(t)];n=ap./(circle.solid ? r : -r)
        else
            ap = [r*cos(t),-r*sin(t)];n = ap./(circle.solid ? r : -r)
        end
    else
        t = asin((f_midpoint[2]-c[2])/r)
        if f_midpoint[1]>c[1]
            ap = [r*cos(t),r*sin(t)];n = ap./ (circle.solid ? r : -r)
        else
            ap = [-r*cos(t),r*sin(t)];n = ap./ (circle.solid ? r : -r)
        end
    end
    return ap,n
end
function calc_intersect(f_midpoint,s_midpoint,::Vector,::Int,circle::Sphere)
    c = circle.center;r = circle.radius
    dir = findfirst(x->x>EPS,abs.(f_midpoint-s_midpoint))
    r_midpoint = f_midpoint-c
    x = r_midpoint[dir%3+1]; y =r_midpoint[(dir+1)%3+1]
    ap = copy(f_midpoint);n = copy(f_midpoint)
    rz = sqrt(r^2-x^2-y^2)
    dz = abs((rz-r_midpoint[dir])/(s_midpoint[dir]-f_midpoint[dir])) < 1 ? rz-r_midpoint[dir] : -rz-r_midpoint[dir]
    ap[dir] = dz+f_midpoint[dir];n[dir] = dz+r_midpoint[dir];n/=(circle.solid ? r : -r)
    return ap,n
end

function pre_partition_box_flag(midpoint,ds,ib::AbstractCircle)
    r = ib.search_radius
    hyper_rec = HyperRectangle(SVector{3,Float64}(ib.center.-ib.radius.-r),SVector{3,Float64}(ib.center.+ib.radius.+r))
    lower = midpoint-0.5*ds;upper = midpoint+0.5*ds
    if overlap_test(lower,upper,hyper_rec)
        return true
    else
        return false
    end
end

function search_radius_refine_flag!(i::Int,ib::AbstractCircle,midpoint,ds,mesh_data)
    r = ib.search_radius
    hyper_rec = HyperRectangle(SVector{3,Float64}(ib.center.-ib.radius.-r),SVector{3,Float64}(ib.center.+ib.radius.+r))
    lower = midpoint-0.5*ds;upper = midpoint+0.5*ds
    if overlap_test(lower,upper,hyper_rec)
        mesh_data.in_box = i
        distance = abs(norm(midpoint-ib.center)-ib.radius)
        if distance<r+0.5*norm(ds)
            mesh_data.in_search_radius=i
            return true
        end
    end
    return false
end

function ghost_cell_flag(ib::Sphere,midpoint,ds)
    flag = 0
    for i = 1:6
        flag += norm(midpoint.+ds.*NMT[3][i].-ib.center)>ib.radius ? 1 : -1 # any neighbor cross boundary?
    end
    abs(flag)==6 && return false
    return true
end