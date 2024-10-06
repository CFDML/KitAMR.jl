function initialize_faces!(
    ::Val{0},
    ::Val{1},
    ip::PW_pxest_iter_face_info_t,
    side,
    amr,
)
    faces = amr.field.faces
    base_quad =
        unsafe_pointer_to_objref(pointer(PointerWrapper(P4est_PS_Data, side.is.full.quad.p.user_data[]).ps_data))
    isa(base_quad,InsideSolidData) && return nothing
    faceid = side.face[] + 1
    any(x->isa(x,AbstractInsideSolidData),base_quad.neighbor.data[faceid]) && return nothing
    push!(faces, Face(InnerFace,base_quad, faceid, nothing))
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
    unsafe_pointer_to_objref(pointer(PointerWrapper(P4est_PS_Data, side.is.full.quad.p.user_data[]).ps_data))
    isa(base_quad,InsideSolidData) && return nothing
    faceid = side.face[] + 1
    any(x->isa(x,AbstractInsideSolidData),base_quad.neighbor.data[faceid]) && return nothing
    push!(faces, Face(InnerFace,base_quad, faceid, nothing))
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
    base_quad =
    unsafe_pointer_to_objref(pointer(PointerWrapper(P4est_PS_Data, qp.p.user_data[]).ps_data))
    isa(base_quad,InsideSolidData) && return nothing
    faceid = side.face[] + 1
    any(x->isa(x,AbstractInsideSolidData),base_quad.neighbor.data[faceid]) && return nothing
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
    push!(faces, Face(InnerFace,base_quad, faceid, hanging_quads))
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
    base_quad =
    unsafe_pointer_to_objref(pointer(PointerWrapper(P4est_PS_Data, qp.p.user_data[]).ps_data))
    isa(base_quad,InsideSolidData) && return nothing
    faceid = side.face[] + 1
    any(x->isa(x,AbstractInsideSolidData),base_quad.neighbor.data[faceid]) && return nothing
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
    push!(faces, Face(InnerFace,base_quad, faceid, hanging_quads))
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
            unsafe_pointer_to_objref(pointer(PointerWrapper(P4est_PS_Data, side.is.full.quad.p.user_data[]).ps_data))
        isa(base_quad,InsideSolidData) && return nothing
        faceid = side.face[] + 1
        any(x->isa(x,AbstractInsideSolidData),base_quad.neighbor.data[faceid]) && return nothing
        push!(faces, Face(InnerFace,base_quad, faceid, nothing))
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
        base_quad =
            unsafe_pointer_to_objref(pointer(PointerWrapper(P4est_PS_Data, qp.p.user_data[]).ps_data))
        isa(base_quad,InsideSolidData) && return nothing
        faceid = side.face[] + 1
        any(x->isa(x,AbstractInsideSolidData),base_quad.neighbor.data[faceid]) && return nothing
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
        push!(faces, Face(InnerFace,base_quad, faceid, hanging_quads))
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
            unsafe_pointer_to_objref(pointer(PointerWrapper(P4est_PS_Data, side.is.full.quad.p.user_data[]).ps_data))
        isa(base_quad,InsideSolidData) && return nothing
        faceid = side.face[] + 1
        any(x->isa(x,AbstractInsideSolidData),base_quad.neighbor.data[faceid]) && return nothing
        push!(faces, Face(InnerFace,base_quad, faceid, nothing))
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
        base_quad =
            unsafe_pointer_to_objref(pointer(PointerWrapper(P4est_PS_Data, qp.p.user_data[]).ps_data))
        isa(base_quad,InsideSolidData) && return nothing
        faceid = side.face[] + 1
        any(x->isa(x,AbstractInsideSolidData),base_quad.neighbor.data[faceid]) && return nothing
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
        push!(faces, Face(InnerFace,base_quad, faceid, hanging_quads))
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
    AMR_face_iterate(info, data, initialize_faces!)
end
function initialize_faces!(ps4est::Ptr{p4est_t}, amr::AMR)
    global_data = amr.global_data
    p_data = pointer_from_objref(amr)
    c_initialize_faces =
        @cfunction(initialize_faces!, Cvoid, (Ptr{p4est_iter_face_info_t}, Ptr{Nothing}))
    GC.@preserve p_data c_initialize_faces AMR_4est_face_iterate(
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
    GC.@preserve p_data c_initialize_faces AMR_4est_face_iterate(
        ps4est,
        global_data.forest.ghost,
        p_data,
        c_initialize_faces,
    )
end

function init_solid_midpoints_kernel(ip, data, dp)
    global_data, solid_midpoints = unsafe_pointer_to_objref(data)
    boundaries = global_data.config.IB
    ds, midpoint = quad_to_cell(ip.p4est, ip.treeid[], ip.quad)
    # flag = true # need to be initialized?
    solid_cell_flags = Vector{Bool}(undef,length(global_data.config.IB))
    for i in eachindex(boundaries)
        inside = solid_flag(boundaries[i],midpoint)
        solid_cell_flags[i] = solid_cell_flag(boundaries[i],midpoint,ds,global_data,inside)
        # flag = flag&&(!inside||solid_cell_flags[i])
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
    global_data, trees, solid_cells = unsafe_pointer_to_objref(data)
    boundaries = global_data.config.IB
    ds, midpoint = quad_to_cell(ip.p4est, ip.treeid[], ip.quad)
    flag = true # need to be initialized?
    solid_cell_flags = Vector{Bool}(undef,length(global_data.config.IB))
    for i in eachindex(boundaries)
        inside = solid_flag(boundaries[i],midpoint)
        solid_cell_flags[i] = solid_cell_flag(boundaries[i],midpoint,ds,global_data,inside)
        flag = flag&&(!inside||solid_cell_flags[i])
        !flag&&break
    end
    treeid = ip.treeid[] - trees.offset
    if flag
        ps_data = PS_Data(typeof(global_data).parameters...)
        push!(trees.data[treeid], ps_data)
        dp[] = P4est_PS_Data(pointer_from_objref(ps_data))
        ic = global_data.config.ic
        gfq = unsafe_wrap(
            Vector{Int},
            pointer(ip.p4est.global_first_quadrant),
            MPI.Comm_size(MPI.COMM_WORLD) + 1,
        )
        ps_data.quadid = global_quadid(ip)+gfq[MPI.Comm_rank(MPI.COMM_WORLD)+1]
        ps_data.ds .= ds
        ps_data.midpoint .= midpoint
        ps_data.prim .= ic
        ps_data.w .= get_conserved(ps_data, global_data)
        ps_data.vs_data = init_VS(ps_data.prim, global_data)
        for i in eachindex(solid_cell_flags)
            if solid_cell_flags[i]
                push!(solid_cells[i].ps_datas,ps_data)
                push!(solid_cells[i].quadids,ps_data.quadid)
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

# function init_PS(comm, global_data)
#     init_PS(comm, global_data, global_data.config.trees_num...)
# end

function re_init_vs4est!(trees, global_data)
    for i in eachindex(trees.data)
        for j in eachindex(trees.data[i])
            ps_data = trees.data[i][j]
            isa(ps_data,InsideSolidData)&&continue
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
    aux_points = init_aux_points(global_data,solid_midpoints)
    broadcast_boundary_midpoints!(aux_points,global_data)
    data = [global_data,aux_points]
    PointerWrapper(ps4est).user_pointer = pointer_from_objref(data)
    GC.@preserve data IB_pre_ps_refine!(ps4est,global_data)
    pre_ps_balance!(ps4est)
    AMR_partition(ps4est)
    return aux_points
end
function IB_Numbers_ranks(ps4est::P_pxest_t,IB_nodes::Vector,solid_cells::Vector{SolidCells{DIM,NDF}}) where{DIM,NDF}
    pp = PointerWrapper(ps4est)
    gfq = unsafe_wrap(
        Vector{Int},
        pointer(pp.global_first_quadrant),
        MPI.Comm_size(MPI.COMM_WORLD) + 1,
    )
    Numbers = Vector{Vector{Vector{Int}}}(undef,MPI.Comm_size(MPI.COMM_WORLD)) # rank{boundaries{solidcells}}
    IB_ranks_table = Vector{Vector{PS_Data{DIM,NDF}}}(undef,MPI.Comm_size(MPI.COMM_WORLD))
    IB_cells = Vector{IBCells}(undef,length(solid_cells))
    rank = MPI.Comm_rank(MPI.COMM_WORLD)
    for i in eachindex(Numbers)
        if i-1 == rank
            Numbers[i] = Vector{Int}[]
        else
            Numbers[i] = Vector{Vector{Int}}(undef,length(solid_cells))
            IB_ranks_table[i] = PS_Data{DIM,NDF}[]
            for j in eachindex(solid_cells)
                Numbers[i][j] = Int[]
                for k in eachindex(solid_cells[j].quadids)
                    if solid_cells[j].quadids[k] >= gfq[i] && solid_cells[j].quadids[k] < gfq[i + 1]
                        push!(Numbers[i][j], length(IB_nodes[j][k])) 
                        append!(IB_ranks_table[i],IB_nodes[j][k])
                    end
                end
            end
        end
    end
    
    for i in eachindex(solid_cells)
        IB_nodes_temp = Vector{Vector{AbstractIBNodes}}(undef,length(solid_cells[i].ps_datas))
        quadids = Vector{Vector{Int}}(undef,length(IB_nodes_temp))
        for j in eachindex(IB_nodes_temp)
            IB_nodes_temp[j] = AbstractIBNodes[]
            quadids[j] = Int[]
        end
        index = 1
        for j in eachindex(solid_cells[i].quadids)
            if solid_cells[i].quadids[j]>= gfq[rank+1]&& solid_cells[i].quadids[j] < gfq[rank + 2]
                append!(IB_nodes_temp[index],IB_nodes[i][j])
                for k in eachindex(IB_nodes[i][j])
                    push!(quadids[index],IB_nodes[i][j][k].quadid)
                end
                index += 1
            end
        end
        IB_cells[i] = IBCells(IB_nodes_temp,quadids)
    end
    return Numbers,IB_ranks_table,IB_cells
end
function IB_nodes_Numbers_communicate(Numbers::Vector{Vector{Vector{Int}}},solid_cells::Vector{T})where{T<:SolidCells}
    rbuffer = Vector{Vector{Int}}(undef,MPI.Comm_size(MPI.COMM_WORLD)) # rank{solid_cells{}} (all boundaries' are together)
    num = 0
    for i in eachindex(solid_cells)
        num+=length(solid_cells[i].ps_datas)
    end
    for i in eachindex(rbuffer)
        i-1==MPI.Comm_rank(MPI.COMM_WORLD)&&continue
        rbuffer[i] = Vector{Int}(undef,num)
    end
    sbuffer = Vector{Vector{Int}}(undef,MPI.Comm_size(MPI.COMM_WORLD))
    for i in eachindex(Numbers)
        i-1 == MPI.Comm_rank(MPI.COMM_WORLD) && continue
        sbuffer[i] = Int[]
        for j in eachindex(Numbers[i])
            append!(sbuffer[i],Numbers[i][j])
        end
    end
    reqs = Vector{MPI.Request}(undef, 0)
    for i in eachindex(sbuffer)
        i-1 == MPI.Comm_rank(MPI.COMM_WORLD)&&continue
        sreq = MPI.Isend(
            sbuffer[i],
            MPI.COMM_WORLD;
            dest = i - 1,
            tag = COMM_DATA_TAG + MPI.Comm_rank(MPI.COMM_WORLD),
        )
        push!(reqs, sreq)
    end
    for i in eachindex(rbuffer)
        i-1 == MPI.Comm_rank(MPI.COMM_WORLD)&&continue
        rreq = MPI.Irecv!(
                rbuffer[i],
                MPI.COMM_WORLD;
                source = i - 1,
                tag = COMM_DATA_TAG + i - 1,
            )
        push!(reqs, rreq)
    end
    MPI.Waitall!(reqs)
    return rbuffer
end
function IB_nodes_vs_nums_communicate(rNumbers::Vector,IB_ranks_table::Vector)
    sbuffer = Vector{Vector{Int}}(undef,MPI.Comm_size(MPI.COMM_WORLD)) # rank{solid_cells{vs_data{}}}
    for i in eachindex(sbuffer)
        i-1 == MPI.Comm_rank(MPI.COMM_WORLD)&&continue
        sbuffer[i] = Vector{Int}[]
        for j in eachindex(IB_ranks_table[i])
            push!(sbuffer[i],IB_ranks_table[i][j].quadid,IB_ranks_table[i][j].vs_data.vs_num)
        end
    end
    reqs = Vector{MPI.Request}(undef, 0)
    for i in eachindex(sbuffer)
        i-1 == MPI.Comm_rank(MPI.COMM_WORLD)&&continue
        isempty(sbuffer[i])&&continue
        sreq = MPI.Isend(
            sbuffer[i],
            MPI.COMM_WORLD;
            dest = i - 1,
            tag = COMM_DATA_TAG + MPI.Comm_rank(MPI.COMM_WORLD),
        )
        push!(reqs, sreq)
    end
    rbuffer = Vector{Vector{Int}}(undef,MPI.Comm_size(MPI.COMM_WORLD)) # ranks{IB_nodes{}}
    r_vs_nums = copy(rbuffer);r_quadids = copy(rbuffer)
    for i in eachindex(rbuffer)
        i-1 == MPI.Comm_rank(MPI.COMM_WORLD)&&continue
        (num = sum(rNumbers[i]))==0 && continue
        rbuffer[i] = Vector{Int}(undef,2*num)
        r_vs_nums[i] = Vector{Int}(undef,num);r_quadids[i] = Vector{Int}(undef,num)
        rreq = MPI.Irecv!(
                rbuffer[i],
                MPI.COMM_WORLD;
                source = i - 1,
                tag = COMM_DATA_TAG + i - 1,
            )
        push!(reqs, rreq)
        MPI.Waitall!(reqs)
        for j in eachindex(r_vs_nums[i])
            r_quadids[i][j] = rbuffer[i][2*j-1]
            r_vs_nums[i][j] = rbuffer[i][2*j]
        end
    end
    return r_quadids,r_vs_nums
end
function IB_nodes_pre_communicate(Numbers::Vector{Vector{Vector{Int}}},solid_cells::Vector{T},IB_ranks_table::Vector) where{T<:SolidCells}
    rNumbers = IB_nodes_Numbers_communicate(Numbers,solid_cells)
    r_quadids,r_vs_nums = IB_nodes_vs_nums_communicate(rNumbers,IB_ranks_table)
    return r_quadids,r_vs_nums
end
function collect_IB_nodes_data(IB_ranks_table::Vector{Vector{PS_Data{DIM,NDF}}}) where{DIM,NDF}
    sdata = Vector{Ptr{Nothing}}(undef,length(IB_ranks_table)) # ranks{solid_cells{}}, data distributes as [df(sum(vs_nums),NDF),level(sum(vs_nums))]
    vnums = Vector{Int}(undef,length(IB_ranks_table))
    IB_buffer = IBBuffer(Ptr{Nothing}[],Ptr{Nothing}[],Vector{Int}[])
    for i in eachindex(IB_ranks_table)
        i-1==MPI.Comm_rank(MPI.COMM_WORLD)&&continue
        isempty(IB_ranks_table[i])&&continue
        for j in eachindex(IB_ranks_table[i])
            vnums[i]+=IB_ranks_table[i][j].vs_data.vs_num
        end
        ps_num = length(IB_ranks_table[i])
        sdata[i] = sc_malloc(-1,vnums[i]*(sizeof(Cdouble)*NDF+sizeof(Int8))+ps_num*DIM*sizeof(Cdouble))
        midpoint_buffer = unsafe_wrap(Matrix{Cdouble},Ptr{Cdouble}(sdata[i]),(ps_num,DIM))
        df_buffer = unsafe_wrap(Matrix{Cdouble},Ptr{Cdouble}(sdata[i]+ps_num*DIM*sizeof(Cdouble)),(vnums[i],NDF))
        level_buffer = unsafe_wrap(Vector{Int8},Ptr{Int8}(sdata[i]+ps_num*DIM*sizeof(Cdouble)+vnums[i]*sizeof(Cdouble)*NDF),vnums[i])
        offset = 0
        for j in eachindex(IB_ranks_table[i])
            midpoint_buffer[j,:] .= IB_ranks_table[i][j].midpoint
        end
        for j in eachindex(IB_ranks_table[i])
            vs_num = IB_ranks_table[i][j].vs_data.vs_num
            df_buffer[offset+1:offset+vs_num,:] .= IB_ranks_table[i][j].vs_data.df
            level_buffer[offset+1:offset+vs_num] .= IB_ranks_table[i][j].vs_data.level
            offset+=vs_num
        end
    end
    IB_buffer.sdata = sdata
    return IB_buffer
end
function IB_nodes_data_communicate(IB_ranks_table::Vector{Vector{PS_Data{DIM,NDF}}},r_vs_nums::Vector) where{DIM,NDF}
    IB_buffer = collect_IB_nodes_data(IB_ranks_table)
    sdata = IB_buffer.sdata
    IB_buffer.rdata = rdata = Vector{Ptr{Nothing}}(undef,MPI.Comm_size(MPI.COMM_WORLD))
    reqs = Vector{MPI.Request}(undef, 0)
    for i in eachindex(sdata)
        i-1 == MPI.Comm_rank(MPI.COMM_WORLD)&&continue
        isempty(IB_ranks_table[i])&&continue
        sreq = MPI.Isend(
            sdata[i],
            MPI.COMM_WORLD;
            dest = i - 1,
            tag = COMM_DATA_TAG + MPI.Comm_rank(MPI.COMM_WORLD),
        )
        push!(reqs, sreq)
    end
    for i in eachindex(rdata)
        i-1 == MPI.Comm_rank(MPI.COMM_WORLD)&&continue
        !isdefined(r_vs_nums,i) && continue
        ps_nums = length(r_vs_nums[i])
        rdata[i] = sc_malloc(-1,sum(r_vs_nums[i])*(sizeof(Cdouble)*NDF+sizeof(Int8))+DIM*ps_nums*sizeof(Cdouble))
        rreq = MPI.Irecv!(
                rdata[i],
                MPI.COMM_WORLD;
                source = i - 1,
                tag = COMM_DATA_TAG + i - 1,
            )
        push!(reqs, rreq)
    end
    MPI.Waitall!(reqs)
    return IB_buffer
end
function IB_nodes_data_wrap!(rNumbers::Vector{Vector{Vector{Int}}},r_quadids::Vector{Vector{Int}},r_vs_nums::Vector{Vector{Int}},rdata::Vector{Ptr{Nothing}},solid_cells::Vector{SolidCells{DIM,NDF}},IB_cells::Vector{IBCells}) where{DIM,NDF}
    rank = MPI.Comm_rank(MPI.COMM_WORLD)
    for i in eachindex(rdata)
        i-1==rank&&continue
        !isdefined(r_vs_nums,i)&&continue
        ps_offset = 1;vs_offset = 0
        vnums = sum(r_vs_nums[i])
        midpoints = unsafe_wrap(Matrix{Cdouble},Ptr{Cdouble}(rdata[i]),(length(r_vs_nums[i]),DIM))
        df_buffer = unsafe_wrap(Matrix{Cdouble},Ptr{Cdouble}(rdata[i]+ps_num*DIM*sizeof(Cdouble)),(vnums,NDF))
        level_buffer = unsafe_wrap(Vector{Int8},Ptr{Int8}(rdata[i]+ps_num*DIM*sizeof(Cdouble)+vnums*sizeof(Cdouble)*NDF),vnums)
        for j in eachindex(IB_cells)
            for k in eachindex(solid_cells[j].ps_datas)
                for _ in 1:rNumbers[i][j][k]
                    vs_num = r_vs_nums[i][ps_offset]
                    push!(IB_cells[j].IB_nodes[k],GhostIBNode{DIM,NDF}(@view(midpoints[ps_offset,:]),@view(level_buffer[vs_offset+1:vs_offset+vs_num]),@view(df_buffer[vs_offset+1:vs_offset+vs_num,:])))
                    push!(IB_cells[j].quadids[k],r_quadids[i][ps_offset])
                    ps_offset+=1;vs_offset+=vs_num
                end
            end
        end
    end
end
function IB_nodes_communicate(Numbers::Vector{Vector{Vector{Int}}},solid_cells::Vector{T},IB_ranks_table::Vector) where{T<:SolidCells}
    r_quadids,r_vs_nums = IB_nodes_pre_communicate(Numbers,solid_cells,IB_ranks_table) # Communicate Numbers and vs_nums for buffer allocation and topological info quadids.
    IB_buffer = IB_nodes_data_communicate(IB_ranks_table,r_vs_nums)
    return r_quadids,r_vs_nums,IB_buffer
end
function init_IBCells(Numbers::Vector{Vector{Vector{Int}}},solid_cells::Vector{T},IB_ranks_table::Vector,IB_cells::Vector{IBCells}) where{T<:SolidCells}
    r_quadids,r_vs_nums,IB_buffer = IB_nodes_communicate(Numbers,solid_cells,IB_ranks_table) # rNumbers: Vector{Vector{Int}}
     IB_nodes_data_wrap!(Numbers,r_quadids,r_vs_nums,IB_buffer.rdata,solid_cells,IB_cells)
     IB_buffer.r_vs_nums = r_vs_nums
     return IB_buffer
end
function init_IB!(ps4est::P_pxest_t,trees::PS_Trees{DIM,NDF},global_data::Global_Data{DIM,NDF},solid_cells::Vector{T},aux_points::Vector) where{DIM,NDF,T<:SolidCells}
    broadcast_quadid!(solid_cells)
    IB_nodes = Vector{Vector{Vector{PS_Data{DIM,NDF}}}}(undef,length(solid_cells)) # boundary{solidcell{}}
    for i in eachindex(IB_nodes)
        IB_nodes[i] = Vector{Vector{PS_Data{DIM,NDF}}}(undef,length(solid_cells[i].quadids))
        for j in eachindex(IB_nodes[i])
            IB_nodes[i][j] = PS_Data{DIM,NDF}[]
        end
    end
    search_IB!(IB_nodes,aux_points,trees,global_data)
    Numbers,IB_ranks_table,IB_cells = IB_Numbers_ranks(ps4est,IB_nodes,solid_cells)
    IB_buffer = init_IBCells(Numbers,solid_cells,IB_ranks_table,IB_cells)
    return IB_cells,IB_buffer,IB_ranks_table
end
function init_ps!(ps4est::P_pxest_t,global_data::Global_Data{DIM,NDF}) where{DIM,NDF}
    fp = PointerWrapper(ps4est)
    trees_data =
        Vector{Vector{AbstractPsData{DIM,NDF}}}(undef, fp.last_local_tree[] - fp.first_local_tree[] + 1)
    for i in eachindex(trees_data)
        trees_data[i] = AbstractPsData{DIM,NDF}[]
    end
    trees = PS_Trees{DIM,NDF}(trees_data, fp.first_local_tree[] - 1)
    solid_cells = Vector{SolidCells{DIM,NDF}}(undef,length(global_data.config.IB))
    for i in eachindex(solid_cells)
        solid_cells[i] = SolidCells{DIM,NDF}(PS_Data{DIM,NDF}[],Int[])
    end
    data = [global_data,trees,solid_cells]
    p_data = pointer_from_objref(data)
    GC.@preserve data AMR_4est_volume_iterate(ps4est, p_data, init_ps_p4est)
    pre_vs_refine!(trees, global_data)
    re_init_vs4est!(trees, global_data)
    return trees, solid_cells
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
        aux_points = pre_refine!(ps4est,global_data)
        trees,solid_cells = init_ps!(ps4est,global_data)
        IB_cells,IB_buffer,IB_ranks_table = init_IB!(ps4est,trees,global_data,solid_cells,aux_points)
        boundary = Boundary{DIM,NDF}(solid_cells,aux_points,IB_cells,IB_ranks_table,IB_buffer)
        return trees, boundary, ps4est
    end
end
function initialize_ghost(p4est::P_pxest_t,global_data::Global_Data)
    ghost_exchange = initialize_ghost_exchange(p4est,global_data)
    ghost_wrap = initialize_ghost_wrap(global_data,ghost_exchange)
    return Ghost(ghost_exchange,ghost_wrap)
end
function init(config::Dict)
    global_data = Global_Data(config)
    trees, boundary, ps4est = init_field!(global_data)
    ghost_ps = AMR_ghost_new(ps4est)
    mesh_ps = AMR_mesh_new(ps4est, ghost_ps)
    global_data.forest.ghost = ghost_ps
    global_data.forest.mesh = mesh_ps
    ghost = initialize_ghost(ps4est, global_data)
    field = Field{config[:DIM],config[:NDF]}(trees,Vector{Face}(undef,0),boundary)
    amr = AMR(
        global_data,ghost,field
    )
    PointerWrapper(ps4est).user_pointer = pointer_from_objref(amr)
    initialize_neighbor_data!(ps4est, amr)
    initialize_faces!(ps4est, amr)
    return (ps4est, amr)
end
