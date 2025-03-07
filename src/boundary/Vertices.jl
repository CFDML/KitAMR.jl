function pre_ps_refine_flag(boundary::Vertices,midpoint::AbstractVector,ds::AbstractVector,global_data::Global_Data) # Circle type IB boundary flag
    boundary.box[1][1]<midpoint[1]&&boundary.box[2][1]>midpoint[1]&&boundary.box[1][2]<midpoint[2]&&boundary.box[2][2]>midpoint[2] && return true
    return false
end
function boundary_flag(boundary::Vertices,midpoint::AbstractVector,ds::AbstractVector,::Global_Data) # Circle type IB boundary flag
    !(boundary.box[1][1]<midpoint[1]&&boundary.box[2][1]>midpoint[1]&&boundary.box[1][2]<midpoint[2]&&boundary.box[2][2]>midpoint[2]) && return false
    flag = 0
    for i = 1:4
        flag += ray_casting(midpoint.+ds.*NMT[2][i],boundary.vertices) ? 1 : -1 # any neighbor cross boundary?
    end
    abs(flag)==4 && return false
    return true
end
function solid_cell_flag(boundary::Vertices,midpoint::AbstractVector,ds::AbstractVector,global_data::Global_Data,inside::Bool) # Ghost nodes, those are inside solid domain and immediately adjacent the boundary.
    (boundary_flag(boundary,midpoint,ds,global_data) && xor(boundary.solid,!inside)) && return true
    return false
end
function solid_flag(boundary::Vertices,midpoint::AbstractVector) # Does midpoint locate at solid?
    inbox = (boundary.box[1][1]<midpoint[1]&&boundary.box[2][1]>midpoint[1]&&boundary.box[1][2]<midpoint[2]&&boundary.box[2][2]>midpoint[2])
    return xor(!(inbox&&ray_casting(midpoint,boundary.vertices)),boundary.solid)
end
function IB_flag(boundary::Vertices,aux_point::AbstractVector,midpoint::AbstractVector,::AbstractVector)
    r = boundary.refine_radius
    norm(midpoint-aux_point)<r
end
function IB_prim(circle::Vertices,aux_point::AbstractVector,ρw::Real)
    IB_prim(circle.bc,aux_point,ρw)
end
function calc_IB_ρw(aux_point::AbstractVector,bound::Vertices{2},midpoint::AbstractMatrix,weight::AbstractVector,df::AbstractMatrix,vn::AbstractVector,Θ::AbstractVector)
    calc_IB_ρw_2D(aux_point,bound.bc,midpoint,weight,df,vn,Θ)
end
function calc_IB_ρw(aux_point::AbstractVector,bound::Vertices{3},midpoint::AbstractMatrix,weight::AbstractVector,df::AbstractMatrix,vn::AbstractVector,Θ::AbstractVector)
    calc_IB_ρw_3D(aux_point,bound.bc,midpoint,weight,df,vn,Θ)
end
function calc_intersect(f_midpoint,s_midpoint,boundary::Vertices)
    find_intersections(s_midpoint,f_midpoint,boundary.vertices)
end