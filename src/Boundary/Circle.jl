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
function _axis_aligned_radial_intersect(f_midpoint, s_midpoint, dir::Int, center, radius, solid::Bool)
    r = Float64(radius)
    radial_sq = r^2
    for d in eachindex(f_midpoint)
        d == dir && continue
        radial_sq -= (f_midpoint[d] - center[d])^2
    end
    if radial_sq < 0.0
        radial_sq < -sqrt(eps(Float64)) * max(r^2, 1.0) &&
            error("No axis-aligned immersed-boundary intersection between fluid center $f_midpoint and solid center $s_midpoint.")
        radial_sq = 0.0
    end

    radial = sqrt(radial_sq)
    candidates = (center[dir] - radial, center[dir] + radial)
    fdir = f_midpoint[dir]
    sdir = s_midpoint[dir]
    dsdir = sdir - fdir
    abs(dsdir) <= EPS && error("Expected an axis-aligned fluid-solid segment in direction $dir.")

    scale = max(abs(fdir), abs(sdir), abs(center[dir]), radial, 1.0)
    tol = sqrt(eps(Float64)) * scale
    lo = min(fdir, sdir) - tol
    hi = max(fdir, sdir) + tol
    best_t = Inf
    best_x = candidates[1]
    for x in candidates
        lo <= x <= hi || continue
        t = clamp((x - fdir) / dsdir, 0.0, 1.0)
        if t < best_t
            best_t = t
            best_x = x
        end
    end
    isfinite(best_t) ||
        error("No immersed-boundary intersection lies between fluid center $f_midpoint and solid center $s_midpoint.")

    ap = copy(f_midpoint)
    ap[dir] = best_x
    n = (ap .- center) ./ (solid ? r : -r)
    return ap, n
end

function calc_intersect(f_midpoint,s_midpoint,::Vector,dir::Int,circle::Circle)
    _axis_aligned_radial_intersect(f_midpoint, s_midpoint, dir, circle.center, circle.radius, circle.solid)
end
function calc_intersect(f_midpoint,s_midpoint,::Vector,dir::Int,circle::Sphere)
    _axis_aligned_radial_intersect(f_midpoint, s_midpoint, dir, circle.center, circle.radius, circle.solid)
end


function search_radius_flag!(i::Int,ib::Circle,midpoint,ds,mesh_data)
    solid_box_flag(midpoint,ds,ib) || return false
    mesh_data.in_box = i
    r = ib.search_radius
    distance = abs(norm(midpoint .- ib.center) - ib.radius)
    if distance < r + 0.5*norm(ds)
        mesh_data.in_search_radius = i
        return true
    end
    return false
end

function search_radius_flag(ib::Circle,midpoint,ds)
    solid_box_flag(midpoint,ds,ib) || return false
    r = ib.search_radius
    distance = abs(norm(midpoint .- ib.center) - ib.radius)
    return distance < r + 0.5*norm(ds)
end

function solid_box_flag(midpoint,ds,ib::Circle)
    r = ib.search_radius
    hyper_rec = HyperRectangle(SVector{2,Float64}(ib.center.-ib.radius.-r),SVector{2,Float64}(ib.center.+ib.radius.+r))
    lower = midpoint-0.5*ds;upper = midpoint+0.5*ds
    if overlap_test(lower,upper,hyper_rec)
        return true
    else
        return false
    end
end
function solid_box_flag(midpoint,ds,ib::Sphere)
    r = ib.search_radius
    hyper_rec = HyperRectangle(SVector{3,Float64}(ib.center.-ib.radius.-r),SVector{3,Float64}(ib.center.+ib.radius.+r))
    lower = midpoint-0.5*ds;upper = midpoint+0.5*ds
    if overlap_test(lower,upper,hyper_rec)
        return true
    else
        return false
    end
end

function search_radius_flag!(i::Int,ib::Sphere,midpoint,ds,mesh_data)
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

"""
$(TYPEDSIGNATURES)
Only accurate for grids that have been refined to the same level (inside the search radius).
"""
function ghost_cell_flag(ib::Circle,midpoint,ds)
    flag = 0
    for i = 1:4
        flag += norm(midpoint.+ds.*NMT[2][i].-ib.center)>ib.radius ? 1 : -1 # any neighbor cross boundary?
    end
    abs(flag)==4 && return false
    return true
end
"""
$(TYPEDSIGNATURES)
Only accurate for grids that have been refined to the same level (inside the search radius).
"""
function ghost_cell_flag(ib::Sphere,midpoint,ds)
    flag = 0
    for i = 1:6
        flag += norm(midpoint.+ds.*NMT[3][i].-ib.center)>ib.radius ? 1 : -1 # any neighbor cross boundary?
    end
    abs(flag)==6 && return false
    return true
end
function search_radius_flag(ib::AbstractCircle,midpoint,ds)
    r = ib.search_radius
    hyper_rec = HyperRectangle(SVector{3,Float64}(ib.center.-ib.radius.-r),SVector{3,Float64}(ib.center.+ib.radius.+r))
    lower = midpoint-0.5*ds;upper = midpoint+0.5*ds
    if overlap_test(lower,upper,hyper_rec)
        distance = abs(norm(midpoint-ib.center)-ib.radius)
        if distance<r+0.5*norm(ds)
            return true
        end
    end
    return false
end
