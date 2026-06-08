# 2D
"""
$(TYPEDSIGNATURES)
"""
function initialize_faces!(p4est::Ptr{p4est_t},ka::KA)
    kinfo = ka.kinfo
    p_data = pointer_from_objref(ka)
    GC.@preserve  ka AMR_face_iterate(p4est;user_data = p_data,ghost = kinfo.forest.ghost) do ip,data
        ka = unsafe_pointer_to_objref(data)
        if ip.sides.elem_count[]==1
            initialize_domain_face!(iPointerWrapper(ip.sides,p4est_iter_face_side_t,0),ka)
        else
            Aside = iPointerWrapper(ip.sides,p4est_iter_face_side_t,0)
            Bside = iPointerWrapper(ip.sides,p4est_iter_face_side_t,1)
            if Aside.is_hanging[]==0
                if Aside.is.full.is_ghost[]==0
                    if Bside.is_hanging[]==0
                        initialize_full_face!(Aside,ka)
                    else
                        initialize_hanging_face!(Aside,ka)
                    end
                elseif Bside.is_hanging[]==0
                    initialize_full_face!(Bside,ka)
                else
                    initialize_back_hanging_face!(Bside,ka)
                end
            elseif Bside.is.full.is_ghost[]==0
                initialize_hanging_face!(Bside,ka)
            else
                initialize_back_hanging_face!(Aside,ka)
            end
        end
    end
end

function solid_full_face_check(base_quad,faceid)
    neighbor = base_quad.neighbor.data[faceid][1]
    (!isa(neighbor,PsData)||neighbor.bound_enc<0) && return true
    return false
end
function solid_hanging_face_check(base_quad,faceid)
    ids = findall(x->(isa(x,PsData)&&x.bound_enc>=0),base_quad.neighbor.data[faceid])
    return ids
end
function initialize_domain_face!(side::PW_pxest_iter_face_side_t,ka::KA{DIM,NDF}) where{DIM,NDF}
    faces = ka.kdata.field.faces
    base_quad = 
        unsafe_pointer_to_objref(pointer(PointerWrapper(P4estPsData,
            side.is.full.quad.p.user_data[]).ps_data))
    isa(base_quad,InsideSolidData) && return nothing
    faceid = side.face[]+1
    direction = get_dir(faceid)
    rot = get_rot(faceid)
    midpoint = copy(base_quad.midpoint)
    midpoint[direction] -= 0.5*rot*base_quad.ds[direction]
    domain = ka.kinfo.config.domain[faceid]
    push!(faces,DomainFace{DIM,NDF,typeof(domain).parameters[1]}(rot,direction,midpoint,domain,base_quad))
    return nothing
end
function initialize_full_face!(side::PW_pxest_iter_face_side_t,ka::KA{DIM,NDF}) where{DIM,NDF}
    faces = ka.kdata.field.faces
    base_quad =
        unsafe_pointer_to_objref(pointer(PointerWrapper(P4estPsData,
            side.is.full.quad.p.user_data[]).ps_data))
    isa(base_quad,InsideSolidData) && return nothing
    faceid = side.face[] + 1
    if base_quad.bound_enc<0
        solid_full_face_check(base_quad,faceid) && return nothing
        base_quad = base_quad.neighbor.data[faceid][1]
		faceid += faceid%2==0 ? -1 : 1
    end
    direction = get_dir(faceid)
    rot = get_rot(faceid)
    midpoint = copy(base_quad.midpoint)
    midpoint[direction] -= 0.5*rot*base_quad.ds[direction]
    kinfo = ka.kinfo
    if midpoint[direction] == kinfo.config.geometry[2*direction-1]||midpoint[direction] == kinfo.config.geometry[2*direction]
        there_midpoint = copy(midpoint);there_midpoint[direction] -= 0.5*rot*base_quad.ds[direction]
        push!(faces,FullFace{DIM,NDF}(rot,direction,midpoint,base_quad,periodic_ghost_cell(there_midpoint,base_quad.neighbor.data[faceid][1])))
    else
        push!(faces,FullFace{DIM,NDF}(rot,direction,midpoint,base_quad,base_quad.neighbor.data[faceid][1]))
    end
    return nothing
end
function initialize_hanging_face!(side::PW_pxest_iter_face_side_t,ka::KA{DIM,NDF}) where{DIM,NDF}
    faces = ka.kdata.field.faces
    base_quad =
        unsafe_pointer_to_objref(pointer(PointerWrapper(P4estPsData,
            side.is.full.quad.p.user_data[]).ps_data))
    isa(base_quad,InsideSolidData) && return nothing
    faceid = side.face[] + 1
    if base_quad.bound_enc<0
        ids = solid_hanging_face_check(base_quad,faceid) 
        isempty(ids) && return nothing     
        here_data = base_quad.neighbor.data[faceid][ids]
        faceid += faceid%2==0 ? -1 : 1
        direction = get_dir(faceid)
        rot = get_rot(faceid)
        midpoint = [copy(x.midpoint) for x in here_data]
        for i in eachindex(midpoint)
            midpoint[i][direction] -= 0.5*rot*here_data[i].ds[direction]
        end
        push!(faces,BackHangingFace{DIM,NDF}(rot,direction,midpoint,here_data,base_quad))
        return nothing
    end
    direction = get_dir(faceid)
    rot = get_rot(faceid)
    neighbor = base_quad.neighbor.data[faceid]
    midpoint = [copy(x.midpoint) for x in neighbor]
    for i in eachindex(neighbor)
        midpoint[i][direction] = base_quad.midpoint[direction] - 0.5 * rot * base_quad.ds[direction]
    end
    kinfo = ka.kinfo
    if  midpoint[1][direction] == kinfo.config.geometry[2*direction-1]||midpoint[1][direction] == kinfo.config.geometry[2*direction]
        periodic_midpoints = [copy(x) for x in midpoint]
        for i in eachindex(periodic_midpoints)
            periodic_midpoints[i][direction] -= 0.5*rot*neighbor[i].ds[direction]
        end
        push!(faces,HangingFace{DIM,NDF}(rot,direction,midpoint,base_quad,periodic_ghost_cell(periodic_midpoints,neighbor)))
    else
        push!(faces,HangingFace{DIM,NDF}(rot,direction,midpoint,base_quad,neighbor))
    end
    return nothing
end
function initialize_back_hanging_face!(side::PointerWrapper{p4est_iter_face_side_t},ka::KA{DIM,NDF}) where{DIM,NDF}
    faces = ka.kdata.field.faces
    is_ghost = Base.unsafe_wrap(
        Vector{Int8},
        Ptr{Int8}(pointer(side.is.hanging.is_ghost)),
        2^(DIM - 1),
    )
    ids = findall(x->x==0,is_ghost).-1
    faceid = side.face[]+1
    here_data = Vector{PsData{DIM,NDF}}()
    for i in ids
        qp = PointerWrapper(iPointerWrapper(side.is.hanging.quad, Ptr{p4est_quadrant_t}, i)[])
        base_quad =
            unsafe_pointer_to_objref(pointer(PointerWrapper(P4estPsData, qp.p.user_data[]).ps_data))
        (!isa(base_quad,PsData)||base_quad.bound_enc<0) && continue
        push!(here_data,base_quad)
    end
    isempty(here_data)&&return nothing
    direction = get_dir(faceid)
    rot = get_rot(faceid)
    midpoint = [copy(x.midpoint) for x in here_data]
    for i in eachindex(midpoint)
        midpoint[i][direction] -= 0.5*rot*here_data[i].ds[direction]
    end
    kinfo = ka.kinfo
    if midpoint[1][direction] == kinfo.config.geometry[2*direction-1]||midpoint[1][direction] == kinfo.config.geometry[2*direction]
        there_midpoint = copy(first(here_data).neighbor.data[faceid][1].midpoint)
        there_midpoint[direction] = midpoint[1][direction] - rot*first(here_data).ds[direction]
        push!(faces,BackHangingFace{DIM,NDF}(rot,direction,midpoint,here_data,periodic_ghost_cell(there_midpoint,first(here_data).neighbor.data[faceid][1])))
    else
        push!(faces,BackHangingFace{DIM,NDF}(rot,direction,midpoint,here_data,first(here_data).neighbor.data[faceid][1]))
    end
    return nothing
end

# 3D
"""
$(TYPEDSIGNATURES)
"""
function initialize_faces!(p4est::Ptr{p8est_t},ka::KA)
    kinfo = ka.kinfo
    p_data = pointer_from_objref(ka)
    GC.@preserve  ka AMR_face_iterate(p4est;user_data = p_data,ghost = kinfo.forest.ghost) do ip,data
        ka = unsafe_pointer_to_objref(data)
        if ip.sides.elem_count[]==1
            initialize_domain_face!(iPointerWrapper(ip.sides,p8est_iter_face_side_t,0),ka)
        else
            Aside = iPointerWrapper(ip.sides,p8est_iter_face_side_t,0)
            Bside = iPointerWrapper(ip.sides,p8est_iter_face_side_t,1)
            if Aside.is_hanging[]==0
                if Aside.is.full.is_ghost[]==0
                    if Bside.is_hanging[]==0
                        initialize_full_face!(Aside,ka)
                    else
                        initialize_hanging_face!(Aside,ka)
                    end
                elseif Bside.is_hanging[]==0
                    initialize_full_face!(Bside,ka)
                else
                    initialize_back_hanging_face!(Bside,ka)
                end
            elseif Bside.is.full.is_ghost[]==0
                initialize_hanging_face!(Bside,ka)
            else
                initialize_back_hanging_face!(Aside,ka)
            end
        end
    end
end
function initialize_back_hanging_face!(side::PointerWrapper{p8est_iter_face_side_t},ka::KA{DIM,NDF}) where{DIM,NDF}
    faces = ka.kdata.field.faces
    is_ghost = Base.unsafe_wrap(
        Vector{Int8},
        Ptr{Int8}(pointer(side.is.hanging.is_ghost)),
        2^(DIM - 1),
    )
    ids = findall(x->x==0,is_ghost).-1
    faceid = side.face[]+1
    here_data = Vector{PsData{DIM,NDF}}()
    for i in ids
        qp = PointerWrapper(iPointerWrapper(side.is.hanging.quad, Ptr{p8est_quadrant_t}, i)[])
        base_quad =
            unsafe_pointer_to_objref(pointer(PointerWrapper(P4estPsData, qp.p.user_data[]).ps_data))
        (!isa(base_quad,PsData)||base_quad.bound_enc<0) && continue
        push!(here_data,base_quad)
    end
    isempty(here_data)&&return nothing
    direction = get_dir(faceid)
    rot = get_rot(faceid)
    midpoint = [copy(x.midpoint) for x in here_data]
        for i in eachindex(midpoint)
            midpoint[i][direction] -= 0.5*rot*here_data[i].ds[direction]
        end
    # push!(faces,BackHangingFace(rot,direction,midpoint,here_data,first(here_data).neighbor.data[faceid][1]))
    kinfo = ka.kinfo
    if midpoint[1][direction] == kinfo.config.geometry[2*direction-1]||midpoint[1][direction] == kinfo.config.geometry[2*direction]
        there_midpoint = copy(first(here_data).neighbor.data[faceid][1].midpoint)
        there_midpoint[direction] = midpoint[1][direction] - rot*first(here_data).ds[direction]
        push!(faces,BackHangingFace{DIM,NDF}(rot,direction,midpoint,here_data,periodic_ghost_cell(there_midpoint,first(here_data).neighbor.data[faceid][1])))
    else
        push!(faces,BackHangingFace{DIM,NDF}(rot,direction,midpoint,here_data,first(here_data).neighbor.data[faceid][1]))
    end
    return nothing
end


function initial_prim(ic::Uniform;kwargs...)
    return ic.ic
end
function initial_prim(ic::PCoordFn;midpoint::AbstractVector{Float64},kinfo::KInfo)
    return ic.PCIC_fn(midpoint,kinfo)
end

function re_init_vs4est!(trees, kinfo)
    for i in eachindex(trees.data)
        for j in eachindex(trees.data[i])
            ps_data = trees.data[i][j]
            isa(ps_data,InsideSolidData)&&continue
            # ps_data.bound_enc<0 && continue
            ps_data.vs_data.df .=
                discrete_maxwell(ps_data, kinfo)
        end
    end
end

"""
$(TYPEDSIGNATURES)
Reapply the configured initial condition on the current physical mesh.

This is intended for initial physical-space AMR: after refinement creates
children from interpolated parent data, call this before rebuilding ghosts so
new fine cells receive the exact initial state at their own cell centers.
"""
function reinitialize_initial_condition!(ka::KA{DIM,NDF}) where{DIM,NDF}
    kinfo = ka.kinfo
    ic = kinfo.config.IC
    for tree in ka.kdata.field.trees.data
        for ps_data in tree
            isa(ps_data, PsData) || continue
            ps_data.bound_enc < 0 && continue
            ps_data.prim .= initial_prim(ic; midpoint = ps_data.midpoint, kinfo = kinfo)
            ps_data.w .= get_conserved(ps_data.prim, kinfo)
            fill!(ps_data.qf, 0.0)
            fill!(ps_data.flux, 0.0)
            fill!(ps_data.sw, 0.0)
            fill!(ps_data.lohner, 0.0)
            ps_data.vs_data = initialize_vs_data(ps_data.prim, kinfo)
        end
    end
    return nothing
end

"""
$(TYPEDSIGNATURES)
"""
function pre_refine!(p4est::Ptr{p4est_t},kinfo::KInfo)
    user_defined_ps_refine!(p4est,kinfo)
    AMR_partition(p4est) do p4est,which_tree,quadrant
        fp = PointerWrapper(p4est)
        kinfo = unsafe_pointer_to_objref(pointer(fp.user_pointer))
        ibs = kinfo.config.IB
        qp = PointerWrapper(quadrant)
        ds,midpoint = quad_to_cell(fp,which_tree,qp)
        for ib in ibs
            if solid_box_flag(midpoint,ds,ib)
                return Cint(1)
            end
        end
        return Cint(0)
    end
    trees = initialize_MeshData!(p4est,kinfo)
    data = [kinfo,trees]
    PointerWrapper(p4est).user_pointer = pointer_from_objref(data)
    GC.@preserve data begin
        search_radius_refine!(p4est,kinfo)
        cell_type_decision!(p4est)
        pre_ps_coarsen!(p4est)
        pre_ps_balance!(p4est)
        meshed_partition!(p4est,trees)
    end
    PointerWrapper(p4est).user_pointer = pointer_from_objref(kinfo)
    return trees
end
"""
$(TYPEDSIGNATURES)
"""
function pre_refine!(p4est::Ptr{p8est_t},kinfo::KInfo)
    user_defined_ps_refine!(p4est,kinfo)
    AMR_partition(p4est) do p4est,which_tree,quadrant
        fp = PointerWrapper(p4est)
        kinfo = unsafe_pointer_to_objref(pointer(fp.user_pointer))
        ibs = kinfo.config.IB
        qp = PointerWrapper(quadrant)
        ds,midpoint = quad_to_cell(fp,which_tree,qp)
        for ib in ibs
            if solid_box_flag(midpoint,ds,ib)
                return Cint(1)
            end
        end
        return Cint(0)
    end
    trees = initialize_MeshData!(p4est,kinfo)
    data = [kinfo,trees]
    PointerWrapper(p4est).user_pointer = pointer_from_objref(data)
    GC.@preserve data begin
        search_radius_refine!(p4est,kinfo)
        cell_type_decision!(p4est)
        pre_ps_coarsen!(p4est)
        pre_ps_balance!(p4est)
        meshed_partition!(p4est,trees)
    end
    PointerWrapper(p4est).user_pointer = pointer_from_objref(kinfo)
    return trees
end

"""
$(TYPEDSIGNATURES)
"""
function initialize_ps!(p4est::Ptr{p4est_t},kinfo::KInfo{DIM,NDF}) where{DIM,NDF}
    fp = PointerWrapper(p4est)
    trees_data =
        Vector{Vector{AbstractPsData{DIM,NDF}}}(undef, fp.last_local_tree[] - fp.first_local_tree[] + 1)
    for i in eachindex(trees_data)
        trees_data[i] = AbstractPsData{DIM,NDF}[]
    end
    trees = PsTrees{DIM,NDF}(trees_data, fp.first_local_tree[] - 1)
    data = [kinfo,trees]
    p_data = pointer_from_objref(data)
    GC.@preserve data AMR_volume_iterate(p4est;user_data = p_data) do ip,data,dp
        kinfo, trees = unsafe_pointer_to_objref(data)
        ds, midpoint = quad_to_cell(ip.p4est, ip.treeid[], ip.quad)
        mesh_data = unsafe_pointer_to_objref(pointer(dp.ps_data))
        treeid = ip.treeid[] - trees.offset # local treeid
        if !mesh_data.in_solid||mesh_data.is_ghost_cell
            ps_data = PsData(typeof(kinfo).parameters...)
            push!(trees.data[treeid], ps_data)
            dp[] = P4estPsData(pointer_from_objref(ps_data))
            ic = kinfo.config.IC
            ps_data.quadid = global_quadid(ip)
            ps_data.ds .= ds
            ps_data.midpoint .= midpoint
            ps_data.prim .= initial_prim(ic;midpoint = ps_data.midpoint,kinfo = kinfo)
            ps_data.w .= get_conserved(ps_data, kinfo)
            ps_data.vs_data = initialize_vs_data(ps_data.prim, kinfo)
            if mesh_data.is_ghost_cell
                ps_data.bound_enc = -mesh_data.in_search_radius
            end
        else
            inside_quad = InsideSolidData{typeof(kinfo).parameters...}(-mesh_data.in_box,midpoint,ds)
            dp[] = P4estPsData(pointer_from_objref(inside_quad))
            push!(trees.data[treeid],inside_quad)
        end
    end
    # re_init_vs4est!(trees, kinfo)
    return trees
end

"""
$(TYPEDSIGNATURES)
"""
function initialize_ps!(p4est::Ptr{p8est_t},kinfo::KInfo{DIM,NDF}) where{DIM,NDF}
    fp = PointerWrapper(p4est)
    trees_data = [AbstractPsData{DIM,NDF}[] for _ in 1:fp.last_local_tree[] - fp.first_local_tree[] + 1]
    trees = PsTrees{DIM,NDF}(trees_data, fp.first_local_tree[] - 1)
    data = [kinfo,trees]
    p_data = pointer_from_objref(data)
    GC.@preserve data AMR_volume_iterate(p4est;user_data = p_data) do ip,data,dp
        kinfo, trees = unsafe_pointer_to_objref(data)
        ds, midpoint = quad_to_cell(ip.p4est, ip.treeid[], ip.quad)
        mesh_data = unsafe_pointer_to_objref(pointer(dp.ps_data))
        treeid = ip.treeid[] - trees.offset # local treeid
        if !mesh_data.in_solid||mesh_data.is_ghost_cell
            ps_data = PsData(typeof(kinfo).parameters...)
            push!(trees.data[treeid], ps_data)
            dp[] = P4estPsData(pointer_from_objref(ps_data))
            ic = kinfo.config.IC
            ps_data.quadid = global_quadid(ip)
            ps_data.ds .= ds
            ps_data.midpoint .= midpoint
            ps_data.prim .= initial_prim(ic;midpoint = ps_data.midpoint,kinfo)
            ps_data.w .= get_conserved(ps_data, kinfo)
            ps_data.vs_data = initialize_vs_data(ps_data.prim, kinfo)
            if mesh_data.is_ghost_cell
                ps_data.bound_enc = -mesh_data.in_search_radius
            end
        else
            inside_quad = InsideSolidData{typeof(kinfo).parameters...}(-mesh_data.in_box,midpoint,ds)
            dp[] = P4estPsData(pointer_from_objref(inside_quad))
            push!(trees.data[treeid],inside_quad)
        end
    end
    # re_init_vs4est!(trees, kinfo)
    return trees
end

"""
$(TYPEDSIGNATURES)
Initialize field for 2D case.
"""
function initialize_trees!(kinfo::KInfo{DIM,NDF}) where{DIM,NDF}
    GC.@preserve kinfo begin
        connectivity_ps = set_connectivity(kinfo)
        p4est = AMR_4est_new(
            MPI.COMM_WORLD,
            connectivity_ps.pointer,
            P4estPsData,
            pointer_from_objref(kinfo),
        )
        kinfo.forest.p4est = p4est
        mesh_tree = pre_refine!(p4est,kinfo)
        GC.@preserve mesh_tree begin
            trees = initialize_ps!(p4est,kinfo)
        end
        return p4est, trees
    end
end

"""
$(TYPEDSIGNATURES)
Initialize field for 3D case.
"""
function initialize_trees!(kinfo::KInfo{3,NDF}) where{NDF}
    GC.@preserve kinfo begin
        connectivity_ps = set_connectivity(kinfo)
        p4est = AMR_4est_new(
            MPI.COMM_WORLD,
            connectivity_ps.pointer,
            P4estPsData,
            pointer_from_objref(kinfo),
        )
        kinfo.forest.p4est = p4est
        mesh_tree = pre_refine!(p4est,kinfo)
        GC.@preserve mesh_tree begin
            trees = initialize_ps!(p4est,kinfo)
        end
        return p4est, trees
    end
end

"""
$(TYPEDSIGNATURES)
Initialize [`Ghost`](@ref) structure.
"""
function initialize_ghost(p4est::P_pxest_t,kinfo::KInfo)
    ghost_buffer, ghost_info = initialize_ghost_pool(p4est,kinfo)
    ghost_wrap = initialize_ghost_wrap(kinfo, ghost_buffer, ghost_info)
    return Ghost(ghost_buffer, ghost_wrap, ghost_info)
end

function initialize_forest!(p4est,kinfo::KInfo)
    kinfo.forest.p4est = p4est
    ghost_ps = AMR_ghost_new(p4est)
    mesh_ps = AMR_mesh_new(p4est, ghost_ps)
    kinfo.forest.ghost = ghost_ps
    kinfo.forest.mesh = mesh_ps
end

"""
$(TYPEDSIGNATURES)
Initial vs_balance. Before calling the function, ghost and neighbor should be initialized.
"""
function initialize_balanced_vs!(ka::KA)
    vs_balance!(ka)
    re_init_vs4est!(ka.kdata.field.trees,ka.kinfo)
    vs_ghost_exchange!(ka.kinfo.forest.p4est, ka)
    update_neighbor!(ka.kinfo.forest.p4est,ka)
end

"""
$(TYPEDSIGNATURES)
Initialize everthing according to `config` dictionary.
"""
function initialize(config::Dict)
    kinfo = KInfo(config)
    p4est, trees = initialize_trees!(kinfo)
    kdata = KData(trees)
    ka = KA(kinfo,kdata)
    PointerWrapper(p4est).user_pointer = pointer_from_objref(ka)
    ps_partition!(p4est, ka)
    initialize_forest!(p4est,kinfo)
    kdata.ghost = initialize_ghost(p4est, kinfo)
    initialize_neighbor_data!(p4est, ka)
    initialize_solid_neighbor!(ka)
    initialize_faces!(p4est, ka)
    initialize_immersed_boundaries!(ka)
    return p4est,ka
end
"""
$(TYPEDSIGNATURES)
Build the full solver state from a [`Configure`](@ref) and return `(p4est, ka)`. This is the
entry point of every run and must be called before any time stepping.

It creates the `p4est` forest, generates and geometry-adaptively refines the physical mesh,
builds the velocity-space grids, initializes the field from the configuration's initial
condition, sets up ghost layers, neighbor maps, faces and immersed boundaries, and performs the
first load-balancing partition.

Returns `(p4est, ka)`:
- `p4est` — opaque pointer to the p4est forest (the parallel mesh topology);
- `ka` — the [`KA`](@ref) object holding all solver data, configuration and status.

Pass both to [`solve!`](@ref) (or to the individual per-step driver functions), then to
[`save_result`](@ref) and [`finalize!`](@ref).
"""
function initialize(config::Configure{DIM,NDF}) where{DIM,NDF}
    kinfo = KInfo(config)
    p4est, trees = initialize_trees!(kinfo)
    kdata = KData(trees)
    ka = KA(kinfo,kdata)
    PointerWrapper(p4est).user_pointer = pointer_from_objref(ka)
    ps_partition!(p4est, ka)
    initialize_forest!(p4est,kinfo)
    kdata.ghost = initialize_ghost(p4est, kinfo)
    initialize_neighbor_data!(p4est, ka)
    initialize_balanced_vs!(ka)
    initialize_solid_neighbor!(ka)
    initialize_faces!(p4est, ka)
    initialize_immersed_boundaries!(ka)
    return p4est,ka
end
