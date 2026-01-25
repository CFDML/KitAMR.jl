function triangle_box_table(kdt::KDTree,mesh::Mesh)
    tree_data = kdt.tree_data
    n_nodes = tree_data.n_internal_nodes
    n_leaves = tree_data.n_leafs
    table = Vector{HyperRectangle{SVector{3,Float64}}}(undef,n_nodes+n_leaves)
    lv = leaves(kdt)
    lower = Inf*ones(3);upper = -Inf*ones(3)
    for node in lv
        id = node.index
        ids = leaf_point_indices(node)
        for j in eachindex(ids)
            vertices = mesh.position[mesh.faces[ids[j]]]
            for k in 1:3
                lower.=min.(lower,vertices[k])
                upper.=max.(upper,vertices[k])
            end
        end
        table[id] = HyperRectangle(SVector{3,Float64}(lower),SVector{3,Float64}(upper))
        lower.=Inf;upper.=-Inf
    end

    for node in PostOrderDFS(treeroot(kdt))
        isempty(children(node))&&continue # is a leaf?
        left_child,right_child = children(node)
        lid = left_child.index;rid = right_child.index
        lower.=min.(table[lid].mins,table[rid].mins)
        upper.=max.(table[lid].maxes,table[rid].maxes)
        table[node.index] = HyperRectangle(SVector{3,Float64}(lower),SVector{3,Float64}(upper))
        lower .= Inf;upper .= -Inf
    end
    return table
end
function triangle_recs(mesh::Mesh)
    recs = Vector{HyperRectangle{SVector{3,Float64}}}(undef,length(mesh.faces))
    lower = Inf*ones(3);upper = -Inf*ones(3)
    for i in eachindex(recs)
        vertices = mesh.position[mesh.faces[i]]
        for j in 1:3
            lower .= min.(lower,vertices[j])
            upper .= max.(upper,vertices[j])
        end
        recs[i] = HyperRectangle(SVector{3,Float64}(lower),SVector{3,Float64}(upper))
        lower .= Inf; upper .= -Inf
    end
    return recs
end
function triangle_edges(mesh::Mesh)
    edges = Vector{Vector{SVector{3,Float64}}}(undef,length(mesh.faces))
    for i in eachindex(mesh.faces)
        vertices = mesh.position[mesh.faces[i]]
        edges[i] = Vector{SVector{3,Float64}}(undef,3)
        for j in 1:3
            edges[i][j] = SVector{3,Float64}(vertices[j%3+1]-vertices[j])
        end
    end
    return edges
end

function overlap_test(lower,upper,hyper_rec::HyperRectangle)
    for i = 1:3
        (upper[i]<hyper_rec.mins[i]||lower[i]>hyper_rec.maxes[i])&& return false
    end
    return true
end
function SAT_test_1(h,C,edge,V0)
    n = cross(edge[1],edge[2])
    l = dot(h,abs.(n))
    s = abs(dot(n,V0-C))
    s>l && return false
    return true
end
function SAT_test_2(a,edge,h,V) # note that V are vertices relative to the midpoint of the cell
    dij = zeros(3);p = zeros(3)
    for j in 1:3
        for i in 1:3
            dij .= cross(a[i],edge[j])
            for k in 1:3
                p[k] = dot(dij,V[k])
            end
            for k in 1:3
                dij[k] = abs(dij[k])
            end
            r = dot(dij,h)
            if minimum(p) > r||maximum(p) < -r # TBD: replace with min/max
                return false
            end
        end
    end
    return true
end



function any_intersect_test(midpoint,ds,tkdt::TriangleKDT)
    lower = midpoint-0.5*ds;upper = midpoint+0.5*ds
    any_recursion_test(treeroot(tkdt.kdt),midpoint,ds,lower,upper,tkdt)
end
function any_recursion_test(node::TreeNode,midpoint,ds,lower,upper,tkdt)
    table = tkdt.table;mesh = tkdt.mesh
    recs = tkdt.triangle_recs;edges = tkdt.triangle_edges
    id = node.index
    if overlap_test(lower,upper,table[id]) # kernel intersect test
        if isempty(children(node))
            ids = leaf_point_indices(node)
            for i in ids
                if overlap_test(lower,upper,recs[i])
                    h = 0.5*ds;V = mesh.position[mesh.faces[i]]
                    if SAT_test_1(h,midpoint,edges[i],V[1])
                        a = [
                            SVector{3,Float64}([1,0,0]),
                            SVector{3,Float64}([0,1,0]),
                            SVector{3,Float64}([0,0,1])
                        ]
                        v = [x-midpoint for x in V]
                        if SAT_test_2(a,edges[i],h,v)
                            return true
                        end
                    end
                end
            end
        else
            lc,rc = children(node)
            any_intersect_recursion_test(lc,midpoint,ds,lower,upper,tkdt) && return true
            return any_intersect_recursion_test(rc,midpoint,ds,lower,upper,tkdt)
        end
    else
        return false
    end
end


function all_intersect_test(midpoint,ds,tkdt::TriangleKDT)
    lower = midpoint-0.5*ds;upper = midpoint+0.5*ds
    intersect_ids = Int[]
    all_intersect_recursion_test!(treeroot(tkdt.kdt),midpoint,ds,lower,upper,tkdt,intersect_ids)
    return intersect_ids
end
function all_intersect_recursion_test!(node::TreeNode,midpoint,ds,lower,upper,tkdt,intersect_ids::Vector{Int})
    table = tkdt.table;mesh = tkdt.mesh
    recs = tkdt.triangle_recs;edges = tkdt.triangle_edges
    id = node.index
    if overlap_test(lower,upper,table[id]) # kernel intersect test
        if isempty(children(node))
            ids = leaf_point_indices(node)
            for i in ids
                if overlap_test(lower,upper,recs[i])
                    h = 0.5*ds;V = mesh.position[mesh.faces[i]]
                    if SAT_test_1(h,midpoint,edges[i],V[1])
                        a = [
                            SVector{3,Float64}([1,0,0]),
                            SVector{3,Float64}([0,1,0]),
                            SVector{3,Float64}([0,0,1])
                        ]
                        v = [x-midpoint for x in V]
                        if SAT_test_2(a,edges[i],h,v)
                            push!(intersect_ids,i)
                        end
                    end
                end
            end
        else
            lc,rc = children(node)
            all_intersect_recursion_test!(lc,midpoint,ds,lower,upper,tkdt,intersect_ids)
            all_intersect_recursion_test!(rc,midpoint,ds,lower,upper,tkdt,intersect_ids)
        end
    end
    return nothing
end



function line_plane_intersection(midpoint::Vector{Float64}, V, edges, dir)
    A,B,C = cross(edges[1],edges[2])
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
function in_triangle_test(point::Vector{Float64},V,edges)
    n = cross(edges[1],edges[2])
    p = [point-v for v in V]
    for i in eachindex(p)
        dot(cross(edges[i],p[i]),n)<-eps()&&return false
    end
    return true
end
function calc_intersect(f_midpoint,s_midpoint,ds,dir,ib::Triangles)
    tkdt = ib.tkdt;mesh = tkdt.mesh
    intersect_ids = all_intersect_test(f_midpoint,2.0*ds,tkdt)
    for id in intersect_ids
        V = mesh.position[mesh.faces[id]]
        edges = ib.tkdt.triangle_edges[id]
        intersect_point = line_plane_intersection(f_midpoint,V,edges,dir)
        isempty(intersect_point) && continue
        (intersect_point[dir]-s_midpoint[dir])*(intersect_point[dir]-f_midpoint[dir])>0 && continue
        abs((intersect_point[dir]-f_midpoint[dir])/ds[dir])>1.0+eps() && continue
        if in_triangle_test(intersect_point,V,edges)
            normal = collect(cross(edges[1],edges[2]))
            normal = dot(normal,f_midpoint-s_midpoint)>0 ? normal/norm(normal) : -normal/norm(normal)
            return intersect_point,normal
        end
    end
    return Float64[],Float64[]
end

function pre_partition_box_flag(midpoint,ds,ib::Triangles)
    r = ib.search_radius
    kdt = ib.tkdt.kdt
    triangle_rec = ib.tkdt.table[treeroot(kdt).index]
    hyper_rec = HyperRectangle(SVector{3,Float64}(triangle_rec.mins.-r),SVector{3,Float64}(triangle_rec.maxes.+r))
    lower = midpoint-0.5*ds;upper = midpoint+0.5*ds
    if overlap_test(lower,upper,hyper_rec)
        return true
    else
        return false
    end
end

function search_radius_refine_flag!(i::Int,ib::Triangles,midpoint,ds,mesh_data)
    r = ib.search_radius
    kdt = ib.tkdt.kdt
    triangle_rec = ib.tkdt.table[treeroot(kdt).index]
    hyper_rec = HyperRectangle(SVector{3,Float64}(triangle_rec.mins.-r),SVector{3,Float64}(triangle_rec.maxes.+r))
    lower = midpoint-0.5*ds;upper = midpoint+0.5*ds
    if overlap_test(lower,upper,hyper_rec)
        mesh_data.in_box = i
        _,distance = nn(kdt,midpoint)
        if distance<r+0.5*norm(ds)
            mesh_data.in_search_radius=i
            return true
        end
    end
    return false
end

function solid_flag(ib::Triangles,midpoint)
    ray_casting(ib,midpoint)
end

function ghost_cell_flag(ib::Triangles,midpoint,ds)
    intersect_ids = all_intersect_test(midpoint,2.0*ds,ib.tkdt)
    isempty(intersect_ids)&&return false
    mesh = ib.tkdt.mesh
    for id in intersect_ids
        V = mesh.position[mesh.faces[id]]
        edges = ib.tkdt.triangle_edges[id]
        for dir in 1:3
            intersect_point = line_plane_intersection(midpoint,V,edges,dir)
            isempty(intersect_point) && continue
            abs((intersect_point[dir]-midpoint[dir])/ds[dir])>1.0+eps() && continue
            in_triangle_test(intersect_point,V,edges) && return true
        end
    end
    return false
end

function cell_type_decision!(p4est::Ptr{p8est_t})
    AMR_volume_iterate(p4est) do ip,data,dp
        global_data,_ = unsafe_pointer_to_objref(pointer(ip.p4est.user_pointer))
        mesh_data = unsafe_pointer_to_objref(pointer(dp.ps_data))
        ds,midpoint = quad_to_cell(ip.p4est,ip.treeid[],ip.quad)
        if mesh_data.in_box!=0
            ib = global_data.config.IB[mesh_data.in_box]
            mesh_data.in_solid = solid_flag(ib,midpoint)
            if mesh_data.in_search_radius!=0&&mesh_data.in_solid
                mesh_data.is_ghost_cell = ghost_cell_flag(ib,midpoint,ds)
            end
        else
            ibs = global_data.config.IB
            if isempty(ibs)
                mesh_data.in_solid = false
            else
                ib = ibs[1]
                mesh_data.in_solid = ib.solid ? false : true
            end
        end
    end
end


function ray_casting(ib::Triangles, midpoint) # Inside the closed surface?
    dir = rand(1:3) # the direction of the ray
    kdt = ib.tkdt.kdt
    center = 0.5*(kdt.hyper_rec.mins[dir]+kdt.hyper_rec.maxes[dir])
    is_pos = center<midpoint[dir]
    intersect_ids = Int[]
    ray_casting_recursion_test!(treeroot(kdt),midpoint,dir,is_pos,ib.tkdt,intersect_ids)
    return isodd(length(intersect_ids))
end
function ray_casting_recursion_test!(node::TreeNode,midpoint,dir::Int,is_pos::Bool,tkdt::TriangleKDT,intersect_ids)
    table = tkdt.table;mesh = tkdt.mesh;recs = tkdt.triangle_recs
    id = node.index
    hyper_rec = table[id]
    if ray_casting_rectangle_test(hyper_rec,midpoint,dir,is_pos)
        if isempty(children(node))
            ids = leaf_point_indices(node)
            for i in ids
                if ray_casting_rectangle_test(recs[i],midpoint,dir,is_pos)
                    V = mesh.position[mesh.faces[i]]
                    if ray_casting_triangle_test(V,midpoint,dir,is_pos)
                        push!(intersect_ids,i)
                    end
                end
            end
        else
            lc,rc = children(node)
            ray_casting_recursion_test!(lc,midpoint,dir,is_pos,tkdt,intersect_ids)
            ray_casting_recursion_test!(rc,midpoint,dir,is_pos,tkdt,intersect_ids)
        end
    end
    return nothing
end
function ray_casting_triangle_test(V,midpoint,dir,is_pos)
    e = zeros(3);e[dir] = is_pos ? 1. : -1.
    Moller_Trumbore_ray_test(midpoint,e,V)
end
function Moller_Trumbore_ray_test(midpoint,e,V)
    e1 = V[2]-V[1];e2 = V[3]-V[1]
    h = cross(e,e2)
    a = dot(h,e1)
    if abs(a)<eps()
        return false
    end
    f = 1.0/a
    s = midpoint-V[1]
    u = f*dot(s,h)
        
    q = cross(s,e1)
    v = f*dot(e,q)

    t = f*dot(e2,q)
    if t>=-eps()&&u>=-eps()&&v>=eps()&&(1.0-u-v>=-eps())
        return true
    else
        return false
    end
end
function Moller_Trumbore_line_test(midpoint,e,V)
    e1 = V[2]-V[1];e2 = V[3]-V[1]
    h = cross(e,e2)
    a = dot(h,e1)
    if abs(a)<eps()
        return false
    end
    f = 1.0/a
    s = midpoint-V[1]
    u = f*dot(s,h)
        
    q = cross(s,e1)
    v = f*dot(e,q)

    if u>=-eps()&&v>=eps()&&(1.0-u-v>=-eps())
        return true
    else
        return false
    end
end
function ray_casting_rectangle_test(hyper_rec::HyperRectangle,midpoint,dir,is_pos)
    dir2 = dir%3+1;dir3 = (dir+1)%3+1
    if (hyper_rec.mins[dir2]<=midpoint[dir2]<=hyper_rec.maxes[dir2])&&(hyper_rec.mins[dir3]<=midpoint[dir3]<=hyper_rec.maxes[dir3])
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

function calc_IB_ρw(aux_point::AbstractVector,bound::Triangles,midpoint::AbstractMatrix,weight::AbstractVector,df::AbstractMatrix,vn::AbstractVector,Θ::AbstractVector)
    calc_IB_ρw_3D(aux_point,bound.bc,midpoint,weight,df,vn,Θ)
end