function pre_ps_refine_flag(boundary::AbstractCircle,midpoint::AbstractVector,ds::AbstractVector,global_data::Global_Data) # Circle type IB boundary flag
    boundary_flag(boundary,midpoint,ds,global_data)
end
function boundary_flag(boundary::Circle,midpoint::AbstractVector,ds::AbstractVector,::Global_Data) # Circle type IB boundary flag
    flag = 0
    for i = 1:4
        flag += norm(midpoint.+ds.*NMT[2][i].-boundary.center)>boundary.radius ? 1 : -1 # any neighbor cross boundary?
    end
    abs(flag)==4 && return false
    return true
end
function boundary_flag(boundary::Sphere,midpoint::AbstractVector,ds::AbstractVector,::Global_Data) # Circle type IB boundary flag
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
function solid_cell_flag(boundary::AbstractCircle,midpoint::AbstractVector,ds::AbstractVector,global_data::Global_Data,inside::Bool) # Ghost nodes, those are inside solid domain and immediately adjacent the boundary.
    (boundary_flag(boundary,midpoint,ds,global_data) && xor(boundary.solid,!inside)) && return true
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
function calc_intersect(f_midpoint,s_midpoint,circle::Circle)
    c = circle.center;r = circle.radius
    if abs(f_midpoint[1]-s_midpoint[1])<EPS
        t = acos((f_midpoint[1]-c[1])/r)
        if f_midpoint[2]>c[2]
            ap = [r*cos(t),r*sin(t)];n=ap./r
        else
            ap = [r*cos(t),-r*sin(t)];n = ap./r
        end
    else
        t = asin((f_midpoint[2]-c[2])/r)
        if f_midpoint[1]>c[1]
            ap = [r*cos(t),r*sin(t)];n = ap./r
        else
            ap = [-r*cos(t),r*sin(t)];n = ap./r
        end
    end
    return ap,n
end