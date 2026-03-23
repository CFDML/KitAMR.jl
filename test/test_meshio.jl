using FileIO, NearestNeighbors, AbstractTrees, StaticArrays, LinearAlgebra
using NearestNeighbors: HyperRectangle
using GeometryBasics: Mesh
mesh = load("./example/X38/X38-surface-coarse.stl")
# face = first(mesh.faces)
centers = zeros(3, length(mesh.faces))
for i in eachindex(mesh.faces)
    vertices = mesh.position[mesh.faces[i]]
    for j = 1:3
        for k = 1:3
            centers[j, i] += vertices[k][j]/3.0
        end
    end
end

kdt = KDTree(centers; leafsize = 25)






# lv = collect(leaves(kdt))
# ids = leaf_point_indices(lv[1])
# # nodes(kdt)
# left_child,right_child = children(treeroot(kdt))
# treeregion(left_child)
# treeregion(right_child)
# treeregion(treeroot(kdt))

function triangle_box_table(kdt::KDTree, mesh::Mesh)
    tree_data = kdt.tree_data
    n_nodes = tree_data.n_internal_nodes
    n_leaves = tree_data.n_leafs
    table = Vector{HyperRectangle{SVector{3,Float64}}}(undef, n_nodes+n_leaves)
    # lv = collect(NearestNeighbors.leaves(kdt))
    lv = leaves(kdt)
    lower = Inf*ones(3);
    upper = -Inf*ones(3)
    for node in lv
        id = node.index
        ids = leaf_point_indices(node)
        for j in eachindex(ids)
            vertices = mesh.position[mesh.faces[ids[j]]]
            for k = 1:3
                lower.=min.(lower, vertices[k])
                upper.=max.(upper, vertices[k])
            end
        end
        table[id] = HyperRectangle(SVector{3,Float64}(lower), SVector{3,Float64}(upper))
        lower.=Inf;
        upper.=-Inf
    end

    for node in PostOrderDFS(treeroot(kdt))
        isempty(children(node))&&continue # is a leaf?
        left_child, right_child = children(node)
        lid = left_child.index;
        rid = right_child.index
        lower.=min.(table[lid].mins, table[rid].mins)
        upper.=max.(table[lid].maxes, table[rid].maxes)
        table[node.index] =
            HyperRectangle(SVector{3,Float64}(lower), SVector{3,Float64}(upper))
        lower .= Inf;
        upper .= -Inf
    end
    return table
end
function triangle_recs(mesh::Mesh)
    recs = Vector{HyperRectangle{SVector{3,Float64}}}(undef, length(mesh.faces))
    lower = Inf*ones(3);
    upper = -Inf*ones(3)
    for i in eachindex(recs)
        vertices = mesh.position[mesh.faces[i]]
        for j = 1:3
            lower .= min.(lower, vertices[j])
            upper .= max.(upper, vertices[j])
        end
        recs[i] = HyperRectangle(SVector{3,Float64}(lower), SVector{3,Float64}(upper))
        lower .= Inf;
        upper .= -Inf
    end
    return recs
end
function triangle_edges(mesh::Mesh)
    edges = Vector{Vector{SVector{3,Float64}}}(undef, length(mesh.faces))
    for i in eachindex(mesh.faces)
        vertices = mesh.position[mesh.faces[i]]
        edges[i] = Vector{SVector{3,Float64}}(undef, 3)
        for j = 1:3
            edges[i][j] = SVector{3,Float64}(vertices[j%3+1]-vertices[j])
        end
    end
    return edges
end

struct TriangleKDT
    kdt::KDTree
    mesh::Mesh
    table::Vector{HyperRectangle{SVector{3,Float64}}}
    triangle_recs::Vector{HyperRectangle{SVector{3,Float64}}}
    triangle_edges::Vector{Vector{SVector{3,Float64}}}
end
function TriangleKDT(mesh::Mesh)
    centers = zeros(3, length(mesh.faces))
    for i in eachindex(mesh.faces)
        vertices = mesh.position[mesh.faces[i]]
        for j = 1:3
            for k = 1:3
                centers[j, i] += vertices[k][j]/3.0
            end
        end
    end
    kdt = KDTree(centers; leafsize = 24)
    return TriangleKDT(
        kdt,
        mesh,
        triangle_box_table(kdt, mesh),
        triangle_recs(mesh),
        triangle_edges(mesh),
    )
end

tkdt = TriangleKDT(mesh)


function any_intersect_test(midpoint, ds, tkdt::TriangleKDT)
    lower = midpoint-0.5*ds;
    upper = midpoint+0.5*ds
    any_recursion_test(treeroot(tkdt.kdt), lower, upper, tkdt)
end
function any_recursion_test(node::TreeNode, lower, upper, tkdt)
    table = tkdt.table;
    mesh = tkdt.mesh
    recs = tkdt.triangle_recs;
    edges = tkdt.triangle_edges
    id = node.index
    if overlap_test(lower, upper, table[id]) # kernel intersect test
        if isempty(children(node))
            ids = leaf_point_indices(node)
            for i in ids
                if overlap_test(lower, upper, recs[i])
                    h = 0.5*ds;
                    V0 = mesh.position[mesh.faces[i]][1]
                    if SAT_test_1(h, midpoint, edges[i], V0)
                        V = [mesh.position[id]-midpoint for id in mesh.faces[i]] # relative position of the vertices
                        a = [
                            SVector{3,Float64}([1, 0, 0]),
                            SVector{3,Float64}([0, 1, 0]),
                            SVector{3,Float64}([0, 0, 1]),
                        ]
                        if SAT_test_2(a, edges[i], h, V)
                            return true
                        end
                    end
                end
            end
        else
            lc, rc = children(node)
            any_intersect_recursion_test(lc, lower, upper, tkdt) && return true
            return any_intersect_recursion_test(rc, lower, upper, tkdt)
        end
    else
        return false
    end
end
function overlap_test(lower, upper, hyper_rec::HyperRectangle)
    for i = 1:3
        (upper[i]<hyper_rec.mins[i]||lower[i]>hyper_rec.maxes[i]) && return false
    end
    return true
end
function SAT_test_1(h, C, edge, V0)
    n = cross(edge[1], edge[2])
    l = dot(h, abs.(n))
    s = abs(dot(n, V0-C))
    s>l && return false
    return true
end
function SAT_test_2(a, edge, h, V) # note that V are vertices relative to the midpoint of the cell
    dij = zeros(3);
    p = zeros(3)
    for j = 1:3
        for i = 1:3
            dij .= cross(a[i], edge[j])
            for k = 1:3
                p[k] = dot(dij, V[k])
            end
            for k = 1:3
                dij[k] = abs(dij[k])
            end
            r = dot(dij, h)
            if minimum(p) > r||maximum(p) < -r # TBD: replace with min/max
                return false
            end
        end
    end
    return true
end



function all_intersect_test(midpoint, ds, tkdt::TriangleKDT)
    lower = midpoint-0.5*ds;
    upper = midpoint+0.5*ds
    intersect_ids = Int[]
    all_intersect_recursion_test!(
        treeroot(tkdt.kdt),
        midpoint,
        ds,
        lower,
        upper,
        tkdt,
        intersect_ids,
    )
    return intersect_ids
end
function all_intersect_recursion_test!(
    node::TreeNode,
    midpoint,
    ds,
    lower,
    upper,
    tkdt,
    intersect_ids::Vector{Int},
)
    table = tkdt.table;
    mesh = tkdt.mesh
    recs = tkdt.triangle_recs;
    edges = tkdt.triangle_edges
    id = node.index
    if overlap_test(lower, upper, table[id]) # kernel intersect test
        if isempty(children(node))
            ids = leaf_point_indices(node)
            for i in ids
                if overlap_test(lower, upper, recs[i])
                    h = 0.5*ds;
                    V = mesh.position[mesh.faces[i]]
                    if SAT_test_1(h, midpoint, edges[i], V[1])
                        a = [
                            SVector{3,Float64}([ds[1], 0, 0]),
                            SVector{3,Float64}([0, ds[2], 0]),
                            SVector{3,Float64}([0, 0, ds[3]]),
                        ]
                        v = [x-midpoint for x in V]
                        if SAT_test_2(a, edges[i], h, v)
                            push!(intersect_ids, i)
                        end
                    end
                end
            end
        else
            lc, rc = children(node)
            all_intersect_recursion_test!(
                lc,
                midpoint,
                ds,
                lower,
                upper,
                tkdt,
                intersect_ids,
            )
            all_intersect_recursion_test!(
                rc,
                midpoint,
                ds,
                lower,
                upper,
                tkdt,
                intersect_ids,
            )
        end
    end
    return nothing
end

midpoint = [5.359375, 0.0390625, 0.90625]
ds = [0.1875, 0.09375, 0.125]
lower = midpoint-0.5ds;
upper = midpoint+0.5ds
intersect_ids = all_intersect_test(midpoint, ds, tkdt)
tri_id = 12209
tri = tkdt.triangle_recs[tri_id]
all(i->tri.mins[i]>lower[i], 1:3)
all(i->tri.maxes[i]<upper[i], 1:3)
mesh = tkdt.mesh
h = 0.5*ds;
V = mesh.position[mesh.faces[tri_id]]
a = [
    SVector{3,Float64}([1, 0, 0]),
    SVector{3,Float64}([0, 1, 0]),
    SVector{3,Float64}([0, 0, 1]),
]
v = [x-midpoint for x in V]
overlap_test(lower, upper, tkdt.triangle_recs[tri_id])
SAT_test_1(h, midpoint, tkdt.triangle_edges[tri_id], V[1])
SAT_test_2(a, tkdt.triangle_edges[tri_id], h, v)

function line_plane_intersection(midpoint::Vector{Float64}, V, edges, dir)
    A, B, C = cross(edges[1], edges[2])
    D = -(A*V[1][1]+B*V[1][2]+C*V[1][3])
    if dir == 1
        y0 = midpoint[2]
        z0 = midpoint[3]
        if A ≈ 0.0
            return Float64[]
        end
        t = -(B*y0 + C*z0 + D) / A
        return [t, y0, z0]
    elseif dir == 2
        x0 = midpoint[1]
        z0 = midpoint[3]
        if B ≈ 0.0
            return Float64[]
        end
        t = -(A*x0 + C*z0 + D) / B
        return [x0, t, z0]
    elseif dir == 3
        x0 = midpoint[1]
        y0 = midpoint[2]
        if C ≈ 0.0
            return Float64[]
        end
        t = -(A*x0 + B*y0 + D) / C
        return [x0, y0, t]
    end
end
function in_triangle_test(point::Vector{Float64}, V, edges)
    n = cross(edges[1], edges[2])
    p = [point-v for v in V]
    for i in eachindex(p)
        dot(cross(edges[i], p[i]), n)<-eps()&&return false
    end
    return true
end

midpoint = [0.25, 0.7999984, 0.0];
ds = ones(3)*0.25
intersect_ids = all_intersect_test(midpoint, ds, tkdt)
for id in intersect_ids
    V = mesh.position[mesh.faces[id]]
    edges = tkdt.triangle_edges[id]
    for dir = 1:3
        intersect_point = line_plane_intersection(midpoint, V, edges, dir)
        isempty(intersect_point) && continue
        abs((intersect_point[dir]-midpoint[dir])/ds[dir])>0.5 && continue
        if in_triangle_test(intersect_point, V, edges)
            @show dir intersect_point
            return true
        end
    end
end





function ray_casting_recursion_test!(
    node::TreeNode,
    midpoint,
    dir::Int,
    is_pos::Bool,
    tkdt::TriangleKDT,
    intersect_ids,
)
    table = tkdt.table;
    mesh = tkdt.mesh;
    edges = tkdt.triangle_edges;
    recs = tkdt.triangle_recs
    id = node.index
    hyper_rec = table[id]
    if ray_casting_rectangle_test(hyper_rec, midpoint, dir, is_pos)
        if isempty(children(node))
            ids = leaf_point_indices(node)
            for i in ids
                if ray_casting_rectangle_test(recs[i], midpoint, dir, is_pos)
                    edge = edges[i]
                    V = mesh.position[mesh.faces[i]]
                    if ray_casting_triangle_test(edge, V, midpoint, dir, is_pos)
                        push!(intersect_ids, i)
                    end
                end
            end
        else
            lc, rc = children(node)
            ray_casting_recursion_test!(lc, midpoint, dir, is_pos, tkdt, intersect_ids)
            ray_casting_recursion_test!(rc, midpoint, dir, is_pos, tkdt, intersect_ids)
        end
    end
    return nothing
end
function ray_casting_triangle_test(edge, V, midpoint, dir, is_pos)
    e = zeros(3);
    e[dir] = is_pos ? 1.0 : -1.0
    h = cross(e, edge[2])
    a = dot(h, edge[1])
    if abs(a)<eps()
        return false
    end
    f = 1.0/a
    s = midpoint-V[1]
    u = f*dot(s, h)
    (u<-eps()||u>1.0+eps())&&return false

    q = cross(s, edge[1])
    v = f*dot(e, q)
    (v<-eps()||v>1.0+eps())&&return false

    t = f*dot(edge[2], q)
    t<eps() && return false

    return true
end
function ray_casting_rectangle_test(hyper_rec::HyperRectangle, midpoint, dir, is_pos)
    dir2 = dir%3+1;
    dir3 = (dir+1)%3+1
    if (
        hyper_rec.mins[dir2]<=midpoint[dir2]<=hyper_rec.maxes[dir2]
    )&&(hyper_rec.mins[dir3]<=midpoint[dir3]<=hyper_rec.maxes[dir3])
        if is_pos
            flag = hyper_rec.maxes[dir]>=midpoint[dir]
        else
            flag = hyper_rec.mins[dir]<=midpoint[dir]
        end
        return flag
    else
        return false
    end
end

function ray_casting(tkdt, midpoint) # Inside the closed surface?
    dir = rand(1:3) # the direction of the ray
    kdt = tkdt.kdt
    center = 0.5*(kdt.hyper_rec.mins[dir]+kdt.hyper_rec.maxes[dir])
    is_pos = center<midpoint[dir]
    intersect_ids = Int[]
    ray_casting_recursion_test!(treeroot(kdt), midpoint, dir, is_pos, tkdt, intersect_ids)
    return isodd(length(intersect_ids))
end
midpoint = [2.0, 1.0, 0.0]
ray_casting(tkdt, midpoint)



mesh = load("./example/X38/X38-surface-fine.stl")
