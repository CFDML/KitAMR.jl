function boundary_flag(
    boundary::Vertices,
    midpoint::AbstractVector,
    ds::AbstractVector,
    ::Union{KInfo,Nothing},
) # Circle type IB boundary flag
    # (boundary.box[1][1]>midpoint[1]||boundary.box[2][1]<midpoint[1]||boundary.box[1][2]>midpoint[2]||boundary.box[2][2]<midpoint[2]) && return false
    flag = 0
    for i = 1:4
        flag += ray_casting(midpoint .+ ds .* NMT[2][i], boundary.vertices) ? 1 : -1 # any neighbor cross boundary?
    end
    abs(flag)==4 && return false
    return true
end
function solid_cell_flag(
    boundary::Vertices,
    midpoint::AbstractVector,
    ds::AbstractVector,
    kinfo::KInfo,
    inside::Bool,
) # Ghost nodes, those are inside solid domain and immediately adjacent the boundary.
    (boundary_flag(boundary, midpoint, ds, kinfo) && inside) && return true
    return false
end
function solid_box_flag(midpoint, ds, ib::Vertices)
    r = ib.search_radius
    lower = midpoint .- 0.5 .* ds
    upper = midpoint .+ 0.5 .* ds
    box_lower = ib.box[1] .- r
    box_upper = ib.box[2] .+ r
    return all(lower .< box_upper) && all(upper .> box_lower)
end
function search_radius_flag!(i::Int, ib::Vertices, midpoint, ds, mesh_data)
    solid_box_flag(midpoint, ds, ib) || return false
    mesh_data.in_box = i
    r = ib.search_radius
    distance = minimum(norm(midpoint .- v) for v in ib.vertices)
    if distance < r + 0.5*norm(ds)
        mesh_data.in_search_radius = i
        return true
    end
    return false
end
function search_radius_flag(ib::Vertices, midpoint, ds)
    solid_box_flag(midpoint, ds, ib) || return false
    r = ib.search_radius
    distance = minimum(norm(midpoint .- v) for v in ib.vertices)
    return distance < r + 0.5*norm(ds)
end
function ghost_cell_flag(ib::Vertices{2}, midpoint, ds)
    return boundary_flag(ib, midpoint, ds, nothing)
end
function solid_flag(boundary::Vertices, midpoint::AbstractVector) # Does midpoint locate at solid?
    inbox = (
        boundary.box[1][1]<midpoint[1]&&boundary.box[2][1]>midpoint[1]&&boundary.box[1][2]<midpoint[2]&&boundary.box[2][2]>midpoint[2]
    )
    return xor(!(inbox&&ray_casting(midpoint, boundary.vertices)), boundary.solid)
end
function IB_prim(circle::Vertices, aux_point::AbstractVector, ρw::Real)
    IB_prim(circle.bc, aux_point, ρw)
end
function calc_IB_ρw(
    aux_point::AbstractVector,
    bound::Vertices{2},
    midpoint::AbstractMatrix,
    weight::AbstractVector,
    df::AbstractMatrix,
    vn::AbstractVector,
    Θ::AbstractVector,
)
    calc_IB_ρw_2D(aux_point, bound.bc, midpoint, weight, df, vn, Θ)
end
function calc_IB_ρw(
    aux_point::AbstractVector,
    bound::Vertices{3},
    midpoint::AbstractMatrix,
    weight::AbstractVector,
    df::AbstractMatrix,
    vn::AbstractVector,
    Θ::AbstractVector,
)
    calc_IB_ρw_3D(aux_point, bound.bc, midpoint, weight, df, vn, Θ)
end
function calc_intersect(f_midpoint, s_midpoint, ::Vector, ::Int, boundary::Vertices{2})
    find_intersections_2D(s_midpoint, f_midpoint, boundary.vertices)
end

function ray_casting(point::Vector{Float64}, vertices::Vector{Vector{Float64}}) # Inside the closed curve?
    n = length(vertices)
    count = 0
    for i = 1:n
        j = (i % n) + 1
        x_i, y_i = vertices[i]
        x_j, y_j = vertices[j]
        px, py = point
        if (py > min(y_i, y_j)) && (py <= max(y_i, y_j))
            if y_i == y_j
                continue
            end
            x_intersect = (py-y_i) * (x_j - x_i) / (y_j - y_i) + x_i
            if x_intersect > px && (x_intersect <= max(x_i, x_j))
                count += 1
            end
        end
    end
    return (count % 2) == 1
end
function calc_normal(p1, p2, s1, s2) # Calculate the unit normal vector of the vector (s2-s1), in the direction of (p2-p1).
    n = [s2[2]-s1[2], s1[1]-s2[1]]
    n /= norm(n)
    if dot(n, p2-p1)<0
        n = -n
    end
    return n
end
function find_horizontal_intersection(s_midpoint, f_midpoint, points)
    n = length(points)
    x_seg_min, x_seg_max = minmax(s_midpoint[1], f_midpoint[1])
    y0 = s_midpoint[2]
    for i = 1:n
        p1 = points[i]
        p2 = points[i%n+1]
        x1, y1 = p1
        x2, y2 = p2
        if (y1 ≤ y0 ≤ y2) || (y2 ≤ y0 ≤ y1)
            t = (y0 - y1) / (y2 - y1)
            if 0 < t < 1
                x_intersect = x1 + t * (x2 - x1)
                if x_seg_min ≤ x_intersect ≤ x_seg_max
                    return [x_intersect, y0], calc_normal(s_midpoint, f_midpoint, p1, p2)
                end
            elseif t==1
                x_intersect = x1 + t * (x2 - x1)
                if x_seg_min ≤ x_intersect ≤ x_seg_max
                    p3 = points[i%n+2]
                    n1 = calc_normal(s_midpoint, f_midpoint, p1, p2)
                    n2 = calc_normal(s_midpoint, f_midpoint, p2, p3)
                    n0 = n1+n2;
                    n0/=norm(n0)
                    return [x_intersect, y0], n0
                end
            end
        end
    end
    @show s_midpoint f_midpoint
    throw("Intersect error!")
    return Float64[], Float64[]
end
function find_vertical_intersection(s_midpoint, f_midpoint, points)
    n = length(points)
    y_seg_min, y_seg_max = minmax(s_midpoint[2], f_midpoint[2])
    x0 = s_midpoint[1]
    for i = 1:n
        p1 = points[i]
        p2 = points[i%n+1]
        x1, y1 = p1
        x2, y2 = p2
        if (x1 ≤ x0 ≤ x2) || (x2 ≤ x0 ≤ x1)
            t = (x0 - x1) / (x2 - x1)
            if 0 < t < 1
                y_intersect = y1 + t * (y2 - y1)
                if y_seg_min ≤ y_intersect ≤ y_seg_max
                    return [x0, y_intersect], calc_normal(s_midpoint, f_midpoint, p1, p2)
                end
            elseif t==1
                y_intersect = y1 + t * (y2 - y1)
                if y_seg_min ≤ y_intersect ≤ y_seg_max
                    p3 = points[i%n+2]
                    n1 = calc_normal(s_midpoint, f_midpoint, p1, p2)
                    n2 = calc_normal(s_midpoint, f_midpoint, p2, p3)
                    n0 = n1+n2;
                    n0/=norm(n0)
                    return [x0, y_intersect], n0
                end
            end
        end
    end
    @show s_midpoint f_midpoint
    throw("Intersect error!")
    return Float64[], Float64[]
end
function find_intersections_2D(s_midpoint, f_midpoint, closed_curve)
    x1, y1 = s_midpoint
    x2, y2 = f_midpoint
    if abs(y1-y2)<EPS
        return find_horizontal_intersection(s_midpoint, f_midpoint, closed_curve)
    elseif abs(x1-x2)<EPS
        return find_vertical_intersection(s_midpoint, f_midpoint, closed_curve)
    else
        @show s_midpoint f_midpoint
        throw("Horizontal or vertical line segment is expected!")
    end
end
