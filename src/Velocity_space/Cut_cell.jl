function gaussian_area(A::AbstractMatrix) # DIMxN
    area = 0
    for i in axes(A,2)
        j = i%size(A,2)+1
        area+=@views det(A[:,[i,j]])
    end
    return 0.5*abs(area)
end

"""
2D axis-aligned cut cell. `gas_weight` is the area in `dot(v,n)<0`; `solid_weight` is
the area in `dot(v,n)>0`.
"""
@inline _positive_square(x::Float64) = x > 0.0 ? x*x : 0.0
function _rect_bounds(vertices::Vector{Vector{Float64}})
    @inbounds begin
        xmin = vertices[1][1];xmax = xmin
        ymin = vertices[1][2];ymax = ymin
        for i in 2:4
            x = vertices[i][1];y = vertices[i][2]
            xmin = min(xmin,x);xmax = max(xmax,x)
            ymin = min(ymin,y);ymax = max(ymax,y)
        end
    end
    return xmin,xmax,ymin,ymax
end
function _cut_rect_negative_halfspace_area(n::Vector{Float64},vertices::Vector{Vector{Float64}})
    xmin,xmax,ymin,ymax = _rect_bounds(vertices)
    Lx = xmax-xmin;Ly = ymax-ymin
    total = Lx*Ly
    total <= 0.0 && return 0.0,total
    n1 = n[1];n2 = n[2]
    ntol = 3.0*eps(Float64)*max(abs(n1),abs(n2))
    a1 = abs(n1) > ntol ? abs(n1) : 0.0
    a2 = abs(n2) > ntol ? abs(n2) : 0.0
    base = 0.0
    a1 > 0.0 && (base += n1 > 0.0 ? n1*xmin : n1*xmax)
    a2 > 0.0 && (base += n2 > 0.0 ? n2*ymin : n2*ymax)
    t = -base
    active = (a1 > 0.0) + (a2 > 0.0)
    active == 0 && return (t >= 0.0 ? total : 0.0),total
    extent = a1*Lx+a2*Ly
    t <= 0.0 && return 0.0,total
    t >= extent && return total,total
    if active == 1
        a1 > 0.0 && return Ly*clamp(t/a1,0.0,Lx),total
        return Lx*clamp(t/a2,0.0,Ly),total
    end
    return (_positive_square(t)-_positive_square(t-a1*Lx)-
        _positive_square(t-a2*Ly)+_positive_square(t-a1*Lx-a2*Ly))/(2.0*a1*a2),total
end
function cut_rect(n::Vector{Float64},vertices::Vector{Vector{Float64}})
    gas_weight,total_weight = _cut_rect_negative_halfspace_area(n,vertices)
    total_weight <= 0.0&&return false,0.,0.
    gas_weight = clamp(gas_weight,0.0,total_weight)
    solid_weight = total_weight-gas_weight
    area_tol = 128.0*eps(Float64)*max(total_weight,1.0)
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
    axis = if abs(n[1]) <= abs(n[2]) && abs(n[1]) <= abs(n[3])
        [1.,0.,0.]
    elseif abs(n[2]) <= abs(n[3])
        [0.,1.,0.]
    else
        [0.,0.,1.]
    end
    e1 = cross(n,axis);e1./=norm(e1)
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
@inline _positive_cube(x::Float64) = x > 0.0 ? x*x*x : 0.0
@inline function _cut_box_volume_2d(t::Float64,a1::Float64,L1::Float64,a2::Float64,L2::Float64,scale::Float64)
    s = _positive_square(t) - _positive_square(t-a1*L1) -
        _positive_square(t-a2*L2) + _positive_square(t-a1*L1-a2*L2)
    return scale*s/(2.0*a1*a2)
end
@inline function _cut_box_volume_3d(t::Float64,a1::Float64,L1::Float64,a2::Float64,L2::Float64,a3::Float64,L3::Float64)
    s = _positive_cube(t) - _positive_cube(t-a1*L1) -
        _positive_cube(t-a2*L2) - _positive_cube(t-a3*L3) +
        _positive_cube(t-a1*L1-a2*L2) + _positive_cube(t-a1*L1-a3*L3) +
        _positive_cube(t-a2*L2-a3*L3) - _positive_cube(t-a1*L1-a2*L2-a3*L3)
    return s/(6.0*a1*a2*a3)
end
function _cut_cube_negative_halfspace_volume(n::Vector{Float64},midpoint::Vector{Float64},ddu::Vector{Float64})
    n1 = n[1];n2 = n[2];n3 = n[3]
    L1 = ddu[1];L2 = ddu[2];L3 = ddu[3]
    total = L1*L2*L3
    ntol = 3.0*eps(Float64)*max(abs(n1),abs(n2),abs(n3))
    a1 = abs(n1) > ntol ? abs(n1) : 0.0
    a2 = abs(n2) > ntol ? abs(n2) : 0.0
    a3 = abs(n3) > ntol ? abs(n3) : 0.0
    base = 0.0
    a1 > 0.0 && (base += n1 > 0.0 ? n1*(midpoint[1]-0.5*L1) : n1*(midpoint[1]+0.5*L1))
    a2 > 0.0 && (base += n2 > 0.0 ? n2*(midpoint[2]-0.5*L2) : n2*(midpoint[2]+0.5*L2))
    a3 > 0.0 && (base += n3 > 0.0 ? n3*(midpoint[3]-0.5*L3) : n3*(midpoint[3]+0.5*L3))
    t = -base
    active = (a1 > 0.0) + (a2 > 0.0) + (a3 > 0.0)
    active == 0 && return t >= 0.0 ? total : 0.0
    extent = a1*L1+a2*L2+a3*L3
    t <= 0.0 && return 0.0
    t >= extent && return total
    if active == 1
        a1 > 0.0 && return L2*L3*clamp(t/a1,0.0,L1)
        a2 > 0.0 && return L1*L3*clamp(t/a2,0.0,L2)
        return L1*L2*clamp(t/a3,0.0,L3)
    elseif active == 2
        a1 == 0.0 && return _cut_box_volume_2d(t,a2,L2,a3,L3,L1)
        a2 == 0.0 && return _cut_box_volume_2d(t,a1,L1,a3,L3,L2)
        return _cut_box_volume_2d(t,a1,L1,a2,L2,L3)
    end
    return _cut_box_volume_3d(t,a1,L1,a2,L2,a3,L3)
end
function cut_cube(n::Vector{Float64},C::Matrix{Float64},midpoint::Vector{Float64},ddu::Vector{Float64},vertices::Matrix{Float64})
    # C and vertices are kept for the existing call signature; the box volume is fixed by midpoint and ddu.
    total_weight = ddu[1]*ddu[2]*ddu[3]
    total_weight <= 0.0 && return false,0.,0.
    gas_weight = clamp(_cut_cube_negative_halfspace_volume(n,midpoint,ddu),0.0,total_weight)
    solid_weight = total_weight-gas_weight
    volume_tol = 128.0*eps(Float64)*max(total_weight,1.0)
    if gas_weight <= volume_tol || solid_weight <= volume_tol
        return false,0.,0.
    end
    return true,gas_weight,solid_weight
end
