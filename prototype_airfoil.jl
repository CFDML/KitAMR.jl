using LinearAlgebra
using CSV
using DataFrames
airfoil  = CSV.read("NACA0012-Original-coordinates.csv", DataFrame; header=true)
airfoil.x
fieldnames(typeof(airfoil))
names(airfoil)

struct Vertices{DIM,T<:AbstractBoundaryType} <: AbstractBoundary
    vertices::Vector{Vector{Float64}} # Vertices of the boundary, sorted in clockwise or counterclockwise order. 
    solid::Bool # Is solid inside the boundary?
    refine_coeffi::Real
    bc::AbstractBCType
    box::Vector{Vector{Float64}} # [[xmin,ymin,zmin],[xmax,ymax,zmax]]
    refine_radius::Real
    Vertices(::Type{T},file::String,solid,search_coeffi,bc) where{T<:AbstractBoundaryType} = (
        source = CSV.read(file,DataFrame;header=true);
        DIM = length(names(source));
        vertices = DIM==2 ? [[s.x[i],s.y[i]] for i in eachindex(source.x)] : [[s.x[i],s.y[i],s.z[i]] for i in eachindex(source.x)];
        box = DIM==2 ? [[minimum(s.x),minimum(s.y)],[maximum(s.x),maximum(s.y)]] : [[minimum(s.x),minimum(s.y),minimum(s.z)],[maximum(s.x),maximum(s.y),maximum(s.z)]];
        new{DIM,T}(vertices,solid,search_coeffi,bc,box)
    )
    Vertices(v::Vertices{DIM,T},config::Configure) where{T<:AbstractBoundaryType}=(
        ds_max = maximum([(config[:geometry][2i]-config[:geometry][2i-1])/config[:trees_num][i] for i in 1:config[:DIM]]);
        ds_min = minimum([(config[:geometry][2i]-config[:geometry][2i-1])/config[:trees_num][i]/2^config[:AMR_PS_MAXLEVEL] for i in 1:config[:DIM]]);
        new{DIM,T}(v.vertices,v.solid,v.search_coeffi,v.bc,[v.box[1].-ds_max,v.box[2].+ds_max],v.search_coeffi*ds_min)
    )
end
# # 待测点坐标
# test_point = (1.5, 0.5)
function ray_casting(point::Vector{Float64},vertices::Vector{Vector{Float64}})
    n = length(vertices)
    count = 0
    for i in 1:n
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
                count +=1
            end
        end
    end
    return (count % 2) == 1
end
vertices = [[airfoil.x[i],airfoil.y[i]] for i in eachindex(airfoil.x)]
ray_casting([0.1,0.02],vertices)
# """
# 判断一点是否在多边形内部（射线法）
# 参数：
# - point: Tuple{Float64, Float64}，待测点的坐标
# - vertices: Vector{Tuple{Float64, Float64}}，多边形的顶点坐标
# 返回值：
# - Bool：true 表示在内部，false 表示在外
# """
# function point_in_polygon(point, vertices)
#     n = length(vertices)
#     count = 0

#     # 遍历每一条边
#     for i in 1:n
#         j = (i % n) + 1  # 下一个顶点的索引

#         # 边的两个端点坐标
#         x_i, y_i = vertices[i]
#         x_j, y_j = vertices[j]

#         # 待测点的坐标
#         px, py = point

#         # 判断射线是否与当前边相交
#         # 射线方向为正x轴，从 (px, py) 沿着 x 增加的方向延伸到无穷远。
#         # 仅当以下条件同时满足时，点和无穷远处之间有一条穿过该边的直线：
#         # 1. 点的 y 坐标在边的两个顶点之间的范围内
#         # 2. 边的斜率不为垂直，且交点在射线的有效区域内

#         # 条件1：点和无穷远处的连线与当前边相交的可能性
#         if (py > min(y_i, y_j)) && (py <= max(y_i, y_j))
#             # 计算交点 x 坐标
#             if x_i == x_j  # 边是垂直于x轴的线段，无法相交（除非点就在该边上）
#                 continue
#             end

#             x_intersect = (y_i - py) * (x_j - x_i) / (y_j - y_i) + x_i

#             # 条件2：交点必须在射线的有效区域内（即 x >= px 且在边的两个顶点之间）
#             if x_intersect > px && (x_intersect <= max(x_i, x_j))
#                 count +=1
#             end
#         end
#     end

#     # 如果交点数为奇数，点在内部；否则在外。
#     return (count % 2) == 1
# end

# # 测试
# result = point_in_polygon(test_point, vertices)
# println("Point $(test_point) is inside the polygon? ", result)
function calc_normal(p1,p2,s1,s2) # Calculate the unit normal vector of the vector (s2-s1), in the direction of (p2-p1).
    n = [s2[2]-s1[2],s1[1]-s2[1]]
    n /= norm(n)
    if dot(n,p2-p1)<0
        n = -n
    end
    return n
end
function find_horizontal_intersection(s_midpoint, f_midpoint, points)
    n = length(points)
    x_seg_min, x_seg_max = minmax(s_midpoint[1], f_midpoint[1])
    y0 = s_midpoint[2]
    for i in 1:n
        p1 = points[i]
        p2 = points[i%n + 1]
        x1, y1 = p1
        x2, y2 = p2
        if (y1 ≤ y0 ≤ y2) || (y2 ≤ y0 ≤ y1)
            t = (y0 - y1) / (y2 - y1)
            if 0 < t < 1
                x_intersect = x1 + t * (x2 - x1)
                if x_seg_min ≤ x_intersect ≤ x_seg_max
                    return [x_intersect, y0],calc_normal(s_midpoint, f_midpoint, p1, p2)
                end
            elseif t==1
                x_intersect = x1 + t * (x2 - x1)
                if x_seg_min ≤ x_intersect ≤ x_seg_max
                    p3 = points[i%n+2]
                    n1 = calc_normal(s_midpoint, f_midpoint, p1, p2)
                    n2 = calc_normal(s_midpoint, f_midpoint, p2, p3)
                    n0 = n1+n2;n0/=norm(n0)
                    return [x_intersect, y0],n0
                end
            end
        end
    end
    throw("Intersect error!")
    return Float64[],Float64[]
end

function find_vertical_intersection(s_midpoint, f_midpoint, points)
    n = length(points)
    y_seg_min, y_seg_max = minmax(s_midpoint[2], f_midpoint[2])
    x0 = s_midpoint[1]
    for i in 1:n
        p1 = points[i]
        p2 = points[i%n + 1]
        x1, y1 = p1
        x2, y2 = p2
        if (x1 ≤ x0 ≤ x2) || (x2 ≤ x0 ≤ x1)
            t = (x0 - x1) / (x2 - x1)
            if 0 < t < 1
                y_intersect = y1 + t * (y2 - y1)
                if y_seg_min ≤ y_intersect ≤ y_seg_max
                    return [x0, y_intersect],calc_normal(s_midpoint, f_midpoint, p1, p2)
                end
            elseif t==1
                y_intersect = y1 + t * (y2 - y1)
                if y_seg_min ≤ y_intersect ≤ y_seg_max
                    p3 = points[i%n+2]
                    n1 = calc_normal(s_midpoint, f_midpoint, p1, p2)
                    n2 = calc_normal(s_midpoint, f_midpoint, p2, p3)
                    n0 = n1+n2;n0/=norm(n0)
                    return [x0, y_intersect],n0
                end
            end
        end
    end
    throw("Intersect error!")
    return Float64[],Float64[]
end

function find_intersections(s_midpoint, f_midpoint, closed_curve)
    x1, y1 = s_midpoint
    x2, y2 = f_midpoint
    if y1 == y2
        y0 = y1
        return find_horizontal_intersection(s_midpoint, f_midpoint, closed_curve)
    elseif x1 == x2
        x0 = x1
        return find_vertical_intersection(s_midpoint, f_midpoint, closed_curve)
    else
        throw("Horizontal or vertical line segment is expected!")
    end
end

source = airfoil
vertices = [[source.x[i],source.y[i]] for i in eachindex(source.x)]
s_midpoint = [0.6720194956400802,-0.036206161383550886];f_midpoint = [0.6720194956400802,-0.03926357056705077]
find_intersections(s_midpoint, f_midpoint, vertices)