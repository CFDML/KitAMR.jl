
function initialize_faces!(ps4est::Ptr{p4est_t},amr::AMR)
    global_data = amr.global_data
    p_data = pointer_from_objref(amr)
    GC.@preserve amr AMR_face_iterate(ps4est;user_data = p_data,ghost = global_data.forest.ghost) do ip,data
        amr = unsafe_pointer_to_objref(data)
        if ip.sides.elem_count[]==1
            initialize_domain_face!(iPointerWrapper(ip.sides,p4est_iter_face_side_t,0),amr)
        else
            Aside = iPointerWrapper(ip.sides,p4est_iter_face_side_t,0)
            Bside = iPointerWrapper(ip.sides,p4est_iter_face_side_t,1)
            if Aside.is_hanging[]==0
                if Aside.is.full.is_ghost[]==0
                    if Bside.is_hanging[]==0
                        initialize_full_face!(Aside,amr)
                    else
                        initialize_hanging_face!(Aside,amr)
                    end
                elseif Bside.is_hanging[]==0
                    initialize_full_face!(Bside,amr)
                else
                    initialize_back_hanging_face!(Bside,amr)
                end
            elseif Bside.is.full.is_ghost[]==0
                initialize_hanging_face!(Bside,amr)
            else
                initialize_back_hanging_face!(Aside,amr)
            end
        end
    end
end
function initialize_domain_face!(side::PointerWrapper{p4est_iter_face_side_t},amr::AMR{DIM,NDF}) where{DIM,NDF}
    faces = amr.field.faces
    base_quad = 
        unsafe_pointer_to_objref(pointer(PointerWrapper(P4est_PS_Data,
            side.is.full.quad.p.user_data[]).ps_data))
    isa(base_quad,InsideSolidData) && return nothing
    faceid = side.face[]+1
    direction = get_dir(faceid)
    rot = get_rot(faceid)
    midpoint = copy(base_quad.midpoint)
    midpoint[direction] -= 0.5*rot*base_quad.ds[direction]
    domain = amr.global_data.config.domain[faceid]
    push!(faces,DomainFace{DIM,NDF,typeof(domain).parameters[1]}(rot,direction,midpoint,domain,base_quad))
    return nothing
end
function solid_full_face_check(base_quad,faceid)
    neighbor = base_quad.neighbor.data[faceid][1]
    (!isa(neighbor,PS_Data)||neighbor.bound_enc<0) && return true
    return false
end
function solid_hanging_face_check(base_quad,faceid)
    ids = findall(x->(isa(x,PS_Data)&&x.bound_enc>=0),base_quad.neighbor.data[faceid])
    return ids
end
function initialize_full_face!(side::PointerWrapper{p4est_iter_face_side_t},amr::AMR{DIM,NDF}) where{DIM,NDF}
    faces = amr.field.faces
    base_quad =
        unsafe_pointer_to_objref(pointer(PointerWrapper(P4est_PS_Data,
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
    push!(faces,FullFace{DIM,NDF}(rot,direction,midpoint,base_quad,base_quad.neighbor.data[faceid][1]))
    return nothing
end
function initialize_hanging_face!(side::PointerWrapper{p4est_iter_face_side_t},amr::AMR{DIM,NDF}) where{DIM,NDF}
    faces = amr.field.faces
    base_quad =
        unsafe_pointer_to_objref(pointer(PointerWrapper(P4est_PS_Data,
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
        midpoint[i][direction] += 0.5*rot*neighbor[i].ds[direction]
    end
    push!(faces,HangingFace{DIM,NDF}(rot,direction,midpoint,base_quad,neighbor))
    return nothing
end
function initialize_back_hanging_face!(side::PointerWrapper{p4est_iter_face_side_t},amr::AMR{DIM,NDF}) where{DIM,NDF}
    faces = amr.field.faces
    is_ghost = Base.unsafe_wrap(
        Vector{Int8},
        Ptr{Int8}(pointer(side.is.hanging.is_ghost)),
        2^(DIM - 1),
    )
    ids = findall(x->x==0,is_ghost).-1
    faceid = side.face[]+1
    here_data = Vector{PS_Data{DIM,NDF}}()
    for i in ids
        qp = PointerWrapper(iPointerWrapper(side.is.hanging.quad, Ptr{p4est_quadrant_t}, i)[])
        base_quad =
            unsafe_pointer_to_objref(pointer(PointerWrapper(P4est_PS_Data, qp.p.user_data[]).ps_data))
        (!isa(base_quad,PS_Data)||base_quad.bound_enc<0) && continue
        push!(here_data,base_quad)
    end
    isempty(here_data)&&return nothing
    direction = get_dir(faceid)
    rot = get_rot(faceid)
    midpoint = [copy(x.midpoint) for x in here_data]
        for i in eachindex(midpoint)
            midpoint[i][direction] -= 0.5*rot*here_data[i].ds[direction]
        end
    push!(faces,BackHangingFace(rot,direction,midpoint,here_data,first(here_data).neighbor.data[faceid][1]))
    return nothing
end

function init_solid_midpoints_kernel(ip, data, dp)
    global_data, solid_midpoints = unsafe_pointer_to_objref(data)
    boundaries = global_data.config.IB
    ds, midpoint = quad_to_cell(ip.p4est, ip.treeid[], ip.quad)
    solid_cell_flags = Vector{Bool}(undef,length(global_data.config.IB))
    for i in eachindex(boundaries)
        inside = solid_flag(boundaries[i],midpoint)
        solid_cell_flags[i] = solid_cell_flag(boundaries[i],midpoint,ds,global_data,inside)
    end
    for i in eachindex(solid_cell_flags)
        if solid_cell_flags[i]
            push!(solid_midpoints[i],midpoint)
        end
    end
end
function init_solid_midpoints(info, data)
    AMR_volume_iterate(info, data, P4est_PS_Data, init_solid_midpoints_kernel)
end

function init_ps_p4est_kernel(ip, data, dp)
    global_data, trees = unsafe_pointer_to_objref(data)
    boundaries = global_data.config.IB
    ds, midpoint = quad_to_cell(ip.p4est, ip.treeid[], ip.quad)
    flag = true # need to be initialized?
    solid_cell_flags = Vector{Bool}(undef,length(global_data.config.IB))
    for i in eachindex(boundaries)
        inside = solid_flag(boundaries[i],midpoint)
	    solid_cell_flags[i] = solid_cell_flag(boundaries[i],midpoint,ds,global_data,inside)&&ip.quad.level[]==global_data.config.solver.AMR_PS_MAXLEVEL
        flag = flag&&(!inside||solid_cell_flags[i])
        !flag&&break
    end
    treeid = ip.treeid[] - trees.offset
    if flag
        ps_data = PS_Data(typeof(global_data).parameters...)
        push!(trees.data[treeid], ps_data)
        dp[] = P4est_PS_Data(pointer_from_objref(ps_data))
        ic = global_data.config.IC
        ps_data.quadid = global_quadid(ip)
        ps_data.ds .= ds
        ps_data.midpoint .= midpoint
        ps_data.prim .= ic
        ps_data.w .= get_conserved(ps_data, global_data)
        ps_data.vs_data = init_VS(ps_data.prim, global_data)
        for i in eachindex(solid_cell_flags)
            if solid_cell_flags[i]
                ps_data.bound_enc<0&&(@error `The solid cell is shared!`)
                ps_data.bound_enc = -i
            end
        end
    else
        inside_quad = InsideSolidData{typeof(global_data).parameters...}()
        dp[] = P4est_PS_Data(pointer_from_objref(inside_quad))
        push!(trees.data[treeid],inside_quad)
    end
end
function init_ps_p4est(info, data)
    AMR_volume_iterate(info, data, P4est_PS_Data, init_ps_p4est_kernel)
end

function re_init_vs4est!(trees, global_data)
    for i in eachindex(trees.data)
        for j in eachindex(trees.data[i])
            ps_data = trees.data[i][j]
            isa(ps_data,InsideSolidData)&&continue
            # ps_data.bound_enc<0 && continue
            ps_data.vs_data.df .=
                discrete_maxwell(ps_data, global_data)
        end
    end
end
function init_aux_points(global_data::Global_Data,solid_midpoints::Vector)
    calc_intersect_point(global_data.config.IB,solid_midpoints)
end
function pre_refine!(ps4est::P_pxest_t,global_data::Global_Data)
    pre_ps_refine!(ps4est,global_data)
    pre_ps_balance!(ps4est)
    solid_midpoints = Vector{Vector{Vector{Float64}}}(undef,length(global_data.config.IB)) # boundaries{solidcells{midpoints{}}}
    for i in eachindex(solid_midpoints)
        solid_midpoints[i] = Vector{Float64}[]
    end
    data = [global_data, solid_midpoints]
    p_data = pointer_from_objref(data)
    GC.@preserve data AMR_4est_volume_iterate(ps4est, p_data, init_solid_midpoints)
    solid_midpoints = broadcast_boundary_midpoints!(solid_midpoints,global_data)
    data = [global_data,solid_midpoints]
    PointerWrapper(ps4est).user_pointer = pointer_from_objref(data)
    GC.@preserve data IB_pre_ps_refine!(ps4est,global_data)
    pre_ps_coarsen!(ps4est;recursive=1)
    pre_ps_balance!(ps4est)
    AMR_partition(ps4est)
end
function init_ps!(ps4est::P_pxest_t,global_data::Global_Data{DIM,NDF}) where{DIM,NDF}
    fp = PointerWrapper(ps4est)
    trees_data =
        Vector{Vector{AbstractPsData{DIM,NDF}}}(undef, fp.last_local_tree[] - fp.first_local_tree[] + 1)
    for i in eachindex(trees_data)
        trees_data[i] = AbstractPsData{DIM,NDF}[]
    end
    trees = PS_Trees{DIM,NDF}(trees_data, fp.first_local_tree[] - 1)
    data = [global_data,trees]
    p_data = pointer_from_objref(data)
    GC.@preserve data AMR_4est_volume_iterate(ps4est, p_data, init_ps_p4est)
    pre_vs_refine!(trees, global_data)
    re_init_vs4est!(trees, global_data)
    return trees
end
function init_field!(global_data::Global_Data{DIM,NDF}) where{DIM,NDF}
    GC.@preserve global_data begin
        connectivity_ps = Cartesian_connectivity(global_data.config.trees_num..., global_data.config.geometry...)
        ps4est = AMR_4est_new(
            MPI.COMM_WORLD,
            connectivity_ps.pointer,
            P4est_PS_Data,
            pointer_from_objref(global_data),
        )
        global_data.forest.p4est = ps4est
        pre_refine!(ps4est,global_data)
        trees = init_ps!(ps4est,global_data)
        return trees, ps4est
    end
end
function initialize_ghost(p4est::P_pxest_t,global_data::Global_Data)
    ghost_exchange = initialize_ghost_exchange(p4est,global_data)
    ghost_wrap = initialize_ghost_wrap(global_data,ghost_exchange)
    return Ghost(ghost_exchange,ghost_wrap)
end
function init(config::Dict)
    global_data = Global_Data(config)
    trees, ps4est = init_field!(global_data)
    ghost_ps = AMR_ghost_new(ps4est)
    mesh_ps = AMR_mesh_new(ps4est, ghost_ps)
    global_data.forest.ghost = ghost_ps
    global_data.forest.mesh = mesh_ps
    ghost = initialize_ghost(ps4est, global_data)
    field = Field{config[:DIM],config[:NDF]}(trees,Vector{AbstractFace}(undef,0))
    amr = AMR(
        global_data,ghost,field
    )
    PointerWrapper(ps4est).user_pointer = pointer_from_objref(amr)
    initialize_neighbor_data!(ps4est, amr)
    initialize_solid_neighbor!(amr)
    reinit_ib_vs!(amr) # Avoid singularity caused by sharp gradient.
    data_exchange!(ps4est, amr)
    initialize_faces!(ps4est, amr)
    return (ps4est, amr)
end
