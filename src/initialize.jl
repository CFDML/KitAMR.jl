function initialize_faces!(
    ::Val{0},
    ::Val{1},
    ip::PW_pxest_iter_face_info_t,
    side,
    amr,
)
    faces = amr.field.faces
    base_quad =
        pointer(PointerWrapper(P4est_PS_Data, side.is.full.quad.p.user_data[]).ps_data)
    faceid = side.face[] + 1
    push!(faces, Face(InnerFace,unsafe_pointer_to_objref(base_quad), faceid, nothing))
end
function initialize_faces!(
    ::Val{0},
    ::Val{0},
    ip::PW_pxest_iter_face_info_t,
    side,
    amr,
)
    faces = amr.field.faces
    base_quad =
        pointer(PointerWrapper(P4est_PS_Data, side.is.full.quad.p.user_data[]).ps_data)
    faceid = side.face[] + 1
    push!(faces, Face(InnerFace,unsafe_pointer_to_objref(base_quad), faceid, nothing))
end
function initialize_faces!(
    ::Val{1},
    ::Val{1},
    ip::PointerWrapper{p4est_iter_face_info_t},
    side::PointerWrapper{T},
    amr::AMR{DIM,NDF},
)   where{T,DIM,NDF}
    faces = amr.field.faces
    side = iPointerWrapper(ip.sides, T, 0)
    is_ghost = Base.unsafe_wrap(
        Vector{Int8},
        Ptr{Int8}(pointer(side.is.hanging.is_ghost)),
        2^(DIM - 1),
    )
    quadid = Base.unsafe_wrap(
        Vector{Int32},
        Ptr{Int32}(pointer(side.is.hanging.quadid)),
        2^(DIM - 1),
    )
    baseid = findfirst(x -> x == 0, is_ghost) - 1
    qp = PointerWrapper(iPointerWrapper(side.is.hanging.quad, Ptr{p4est_quadrant_t}, baseid)[])
    base_quad = pointer(PointerWrapper(P4est_PS_Data, qp.p.user_data[]).ps_data)
    faceid = side.face[] + 1
    index = 1
    hanging_quads = Vector{AbstractPsData{DIM,NDF}}(undef, 2^(DIM - 1)-1)
    for i in 0:2^(DIM-1)-1
        i==baseid && continue
        qp = PointerWrapper(
            iPointerWrapper(side.is.hanging.quad, Ptr{p4est_quadrant_t}, i)[],
        )
        if is_ghost[i+1] == 0
            hanging_quads[index] = unsafe_pointer_to_objref(
                pointer(PointerWrapper(P4est_PS_Data, qp.p.user_data[]).ps_data),
            )
        elseif quadid[i+1]>-1
            hanging_quads[index] = amr.ghost.ghost_wrap[quadid[i+1]+1]
        else
            hanging_quads[index] = MissingHangingQuad{DIM,NDF}()
        end
        index+=1
    end
    push!(faces, Face(InnerFace,unsafe_pointer_to_objref(base_quad), faceid, hanging_quads))
end
function initialize_faces!(
    ::Val{1},
    ::Val{1},
    ip::PointerWrapper{p8est_iter_face_info_t},
    side::PointerWrapper{T},
    amr::AMR{DIM,NDF},
)   where{T,DIM,NDF}
    faces = amr.field.faces
    side = iPointerWrapper(ip.sides, T, 0)
    is_ghost = Base.unsafe_wrap(
        Vector{Int8},
        Ptr{Int8}(pointer(side.is.hanging.is_ghost)),
        2^(DIM - 1),
    )
    quadid = Base.unsafe_wrap(
        Vector{Int32},
        Ptr{Int32}(pointer(side.is.hanging.quadid)),
        2^(DIM - 1),
    )
    baseid = findfirst(x -> x == 0, is_ghost) - 1
    qp = PointerWrapper(iPointerWrapper(side.is.hanging.quad, Ptr{p8est_quadrant_t}, baseid)[])
    base_quad = pointer(PointerWrapper(P4est_PS_Data, qp.p.user_data[]).ps_data)
    faceid = side.face[] + 1
    index = 1
    hanging_quads = Vector{AbstractPsData{DIM,NDF}}(undef, 2^(DIM - 1)-1)
    for i in 0:2^(DIM-1)-1
        i==baseid && continue
        qp = PointerWrapper(
            iPointerWrapper(side.is.hanging.quad, Ptr{p8est_quadrant_t}, i)[],
        )
        if is_ghost[i+1] == 0
            hanging_quads[index] = unsafe_pointer_to_objref(
                pointer(PointerWrapper(P4est_PS_Data, qp.p.user_data[]).ps_data),
            )
        elseif quadid[i+1]>-1
            hanging_quads[index] = amr.ghost.ghost_wrap[quadid[i+1]+1]
        else
            hanging_quads[index] = MissingHangingQuad{DIM,NDF}()
        end
        index+=1
    end
    push!(faces, Face(InnerFace,unsafe_pointer_to_objref(base_quad), faceid, hanging_quads))
end
function initialize_faces!(
    ::Val{1},
    ::Val{0},
    ip::PointerWrapper{p4est_iter_face_info_t},
    side::PointerWrapper{T},
    amr::AMR{DIM,NDF},
) where{T,DIM,NDF}
    side = iPointerWrapper(ip.sides, T, 1)
    faces = amr.field.faces
    if side.is_hanging[] == 0
        base_quad =
            pointer(PointerWrapper(P4est_PS_Data, side.is.full.quad.p.user_data[]).ps_data)
        faceid = side.face[] + 1
        push!(faces, Face(unsafe_pointer_to_objref(base_quad), faceid, 0, nothing))
    else
        is_ghost = Base.unsafe_wrap(
            Vector{Int8},
            Ptr{Int8}(pointer(side.is.hanging.is_ghost)),
            2^(DIM - 1),
        )
        quadid = Base.unsafe_wrap(
            Vector{Int32},
            Ptr{Int32}(pointer(side.is.hanging.quadid)),
            2^(DIM - 1),
        )
        baseid = findfirst(x -> x == 0, is_ghost) - 1
        qp = PointerWrapper(
            iPointerWrapper(side.is.hanging.quad, Ptr{p4est_quadrant_t}, baseid)[],
        )
        base_quad = pointer(PointerWrapper(P4est_PS_Data, qp.p.user_data[]).ps_data)
        faceid = side.face[] + 1
        index = 1
        hanging_quads = Vector{AbstractPsData{DIM,NDF}}(undef, 2^(DIM - 1)-1)
        for i in 0:2^(DIM-1)-1
            i==baseid && continue
            qp = PointerWrapper(
                iPointerWrapper(side.is.hanging.quad, Ptr{p4est_quadrant_t}, i)[],
            )
            if is_ghost[i+1] == 0
                hanging_quads[index] = unsafe_pointer_to_objref(
                    pointer(PointerWrapper(P4est_PS_Data, qp.p.user_data[]).ps_data),
                )
            elseif quadid[i+1]>-1
                hanging_quads[index] = amr.ghost.ghost_wrap[quadid[i+1]+1]
            else
                hanging_quads[index] = MissingHangingQuad{DIM,NDF}()
            end
            index+=1
        end                
        push!(faces, Face(InnerFace,unsafe_pointer_to_objref(base_quad), faceid, hanging_quads))
    end
end
function initialize_faces!(
    ::Val{1},
    ::Val{0},
    ip::PointerWrapper{p8est_iter_face_info_t},
    side::PointerWrapper{T},
    amr::AMR{DIM,NDF},
) where{T,DIM,NDF}
    side = iPointerWrapper(ip.sides, T, 1)
    faces = amr.field.faces
    if side.is_hanging[] == 0
        base_quad =
            pointer(PointerWrapper(P4est_PS_Data, side.is.full.quad.p.user_data[]).ps_data)
        faceid = side.face[] + 1
        push!(faces, Face(InnerFace,unsafe_pointer_to_objref(base_quad), faceid, nothing))
    else
        is_ghost = Base.unsafe_wrap(
            Vector{Int8},
            Ptr{Int8}(pointer(side.is.hanging.is_ghost)),
            2^(DIM - 1),
        )
        quadid = Base.unsafe_wrap(
            Vector{Int32},
            Ptr{Int32}(pointer(side.is.hanging.quadid)),
            2^(DIM - 1),
        )
        baseid = findfirst(x -> x == 0, is_ghost) - 1
        qp = PointerWrapper(
            iPointerWrapper(side.is.hanging.quad, Ptr{p8est_quadrant_t}, baseid)[],
        )
        base_quad = pointer(PointerWrapper(P4est_PS_Data, qp.p.user_data[]).ps_data)
        faceid = side.face[] + 1
        index = 1
        hanging_quads = Vector{AbstractPsData{DIM,NDF}}(undef, 2^(DIM - 1)-1)
        for i in 0:2^(DIM-1)-1
            i==baseid && continue
            qp = PointerWrapper(
                iPointerWrapper(side.is.hanging.quad, Ptr{p8est_quadrant_t}, i)[],
            )
            if is_ghost[i+1] == 0
                hanging_quads[index] = unsafe_pointer_to_objref(
                    pointer(PointerWrapper(P4est_PS_Data, qp.p.user_data[]).ps_data),
                )
            elseif quadid[i+1]>-1
                hanging_quads[index] = amr.ghost.ghost_wrap[quadid[i+1]+1]
            else
                hanging_quads[index] = MissingHangingQuad{DIM,NDF}()
            end
            index+=1
        end                
        push!(faces, Face(InnerFace,unsafe_pointer_to_objref(base_quad), faceid, hanging_quads))
    end
end
function initialize_faces!(
    ::Val{2},
    ip::PointerWrapper{p4est_iter_face_info_t},
    data::Ptr{Nothing},
)
    amr = unsafe_pointer_to_objref(data)
    sideid = 0
    side = iPointerWrapper(ip.sides, p4est_iter_face_side_t, sideid)
    side.is_hanging[] == 1 &&
        (side = iPointerWrapper(ip.sides, p4est_iter_face_side_t, (sideid += 1)))
    initialize_faces!(Val(Int(side.is.full.is_ghost[])), Val(sideid), ip, side, amr)
end
function initialize_faces!(
    ::Val{2},
    ip::PointerWrapper{p8est_iter_face_info_t},
    data::Ptr{Nothing},
)
    amr = unsafe_pointer_to_objref(data)
    sideid = 0
    side = iPointerWrapper(ip.sides, p8est_iter_face_side_t, sideid)
    side.is_hanging[] == 1 &&
        (side = iPointerWrapper(ip.sides, p8est_iter_face_side_t, (sideid += 1)))
    initialize_faces!(Val(Int(side.is.full.is_ghost[])), Val(sideid), ip, side, amr)
end
function initialize_faces!(
    ::Val{1},
    ip::PointerWrapper{p4est_iter_face_info_t},
    data::Ptr{Nothing},
)
    amr = unsafe_pointer_to_objref(data)
    faces = amr.field.faces
    side = iPointerWrapper(ip.sides, p4est_iter_face_side_t, 0)
    base_quad =
        pointer(PointerWrapper(P4est_PS_Data, side.is.full.quad.p.user_data[]).ps_data)
    faceid = side.face[] + 1
    push!(faces, Face(BoundaryFace,unsafe_pointer_to_objref(base_quad), faceid, nothing))
end
function initialize_faces!(
    ::Val{1},
    ip::PointerWrapper{p8est_iter_face_info_t},
    data::Ptr{Nothing},
)
    amr = unsafe_pointer_to_objref(data)
    faces = amr.field.faces
    side = iPointerWrapper(ip.sides, p8est_iter_face_side_t, 0)
    base_quad =
        pointer(PointerWrapper(P4est_PS_Data, side.is.full.quad.p.user_data[]).ps_data)
    faceid = side.face[] + 1
    push!(faces, Face(BoundaryFace,unsafe_pointer_to_objref(base_quad), faceid, nothing))
end
function initialize_faces!(ip::PW_pxest_iter_face_info_t, data)
    initialize_faces!(Val(Int(ip.sides.elem_count[])), ip, data)
end
function initialize_faces!(info::P_pxest_iter_face_info_t, data)
    DVM_face_iterate(info, data, initialize_faces!)
end
function initialize_faces!(ps4est::Ptr{p4est_t}, amr::AMR)
    global_data = amr.global_data
    p_data = pointer_from_objref(amr)
    c_initialize_faces =
        @cfunction(initialize_faces!, Cvoid, (Ptr{p4est_iter_face_info_t}, Ptr{Nothing}))
    GC.@preserve p_data c_initialize_faces DVM_4est_face_iterate(
        ps4est,
        global_data.forest.ghost,
        p_data,
        c_initialize_faces,
    )
end
function initialize_faces!(ps4est::Ptr{p8est_t}, amr::AMR)
    global_data = amr.global_data
    p_data = pointer_from_objref(amr)
    c_initialize_faces =
        @cfunction(initialize_faces!, Cvoid, (Ptr{p8est_iter_face_info_t}, Ptr{Nothing}))
    GC.@preserve p_data c_initialize_faces DVM_4est_face_iterate(
        ps4est,
        global_data.forest.ghost,
        p_data,
        c_initialize_faces,
    )
end

function init_p4est_kernal(ip, data, dp)
    global_data, trees = unsafe_pointer_to_objref(data)
    treeid = ip.treeid[] - trees.offset
    ps_data = PS_Data(typeof(global_data).parameters...)
    push!(trees.data[treeid], ps_data)
    dp[] = P4est_PS_Data(pointer_from_objref(ps_data))
    ic = global_data.config.ic
    ds, midpoint = quad_to_cell(ip.p4est, ip.treeid[], ip.quad)
    ps_data.ds .= ds
    ps_data.midpoint .= midpoint
    ps_data.prim .= ic
    ps_data.w .= get_conserved(ps_data, global_data)
    ps_data.vs_data = init_VS(ps_data.prim, global_data)
end
function init_ps4est(info, data)
    DVM_volume_iterate(info, data, P4est_PS_Data, init_p4est_kernal)
end
function init_PS(comm, global_data)
    init_PS(comm, global_data, global_data.config.trees_num...)
end
function re_init_vs4est!(trees, global_data)
    for i in eachindex(trees.data)
        for j in eachindex(trees.data[i])
            ps_data = trees.data[i][j]
            ps_data.vs_data.df .=
                discrete_maxwell(ps_data, global_data)
        end
    end
end
function init_PS(comm::MPI.Comm, global_data::Global_Data{DIM,NDF}) where{DIM,NDF}
    GC.@preserve global_data begin
        connectivity_ps = Cartesian_connectivity(global_data.config.trees_num..., global_data.config.geometry...)
        ps4est = DVM_4est_new(
            comm,
            connectivity_ps.pointer,
            P4est_PS_Data,
            pointer_from_objref(global_data),
        )
        pre_ps_refine!(ps4est,global_data)
        pre_ps_balance!(ps4est)
        AMR_partition(ps4est)
        fp = PointerWrapper(ps4est)
        trees_data =
            Vector{Vector{PS_Data{DIM,NDF}}}(undef, fp.last_local_tree[] - fp.first_local_tree[] + 1)
        for i in eachindex(trees_data)
            trees_data[i] = Vector{PS_Data{DIM,NDF}}(undef, 0)
        end
        trees = PS_Trees{DIM,NDF}(trees_data, fp.first_local_tree[] - 1)
        data = [global_data, trees]
        p_data = pointer_from_objref(data)
        GC.@preserve data DVM_4est_volume_iterate(ps4est, p_data, init_ps4est)
        pre_vs_refine!(trees, global_data)
        re_init_vs4est!(trees, global_data)
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
    trees, ps4est = init_PS(MPI.COMM_WORLD, global_data)
    ghost_ps = AMR_ghost_new(ps4est)
    mesh_ps = AMR_mesh_new(ps4est, ghost_ps)
    global_data.forest.ghost = ghost_ps
    global_data.forest.mesh = mesh_ps
    ghost = initialize_ghost(ps4est, global_data)
    field = Field{config[:DIM],config[:NDF]}(trees,Vector{Face}(undef,0))
    amr = AMR(
        global_data,ghost,field
    )
    initialize_faces!(ps4est, amr)
    PointerWrapper(ps4est).user_pointer = pointer_from_objref(amr)
    initialize_neighbor_data!(ps4est, amr)
    return (ps4est, amr)
end
