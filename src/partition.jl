function partition_check(p4est::P_pxest_t)
    pp = PointerWrapper(p4est)
    gfq = Base.unsafe_wrap(
        Vector{Int},
        pointer(pp.global_first_quadrant),
        MPI.Comm_size(MPI.COMM_WORLD) + 1,
    )
    nums = [gfq[i]-gfq[i-1] for i in 2:MPI.Comm_size(MPI.COMM_WORLD)+1]
    return (maximum(nums)-minimum(nums))/minimum(nums)>0.1
end
function partition!(p4est::Ptr{p4est_t})
    pp = PointerWrapper(p4est)
    gfq = Base.unsafe_wrap(
        Vector{Int},
        pointer(pp.global_first_quadrant),
        MPI.Comm_size(MPI.COMM_WORLD) + 1,
    )
    src_gfq = copy(gfq)
    src_flt = pp.first_local_tree[]
    src_llt = pp.last_local_tree[]
    p4est_partition(p4est, 0, C_NULL)
    # dest_gfq = Base.unsafe_wrap(Vector{Int},pointer(pp.global_first_quadrant),MPI.Comm_size(MPI.COMM_WORLD)+1)
    return src_gfq, gfq, src_flt, src_llt
end
function partition!(p4est::Ptr{p8est_t})
    pp = PointerWrapper(p4est)
    gfq = Base.unsafe_wrap(
        Vector{Int},
        pointer(pp.global_first_quadrant),
        MPI.Comm_size(MPI.COMM_WORLD) + 1,
    )
    src_gfq = copy(gfq)
    src_flt = pp.first_local_tree[]
    src_llt = pp.last_local_tree[]
    p8est_partition(p4est, 0, C_NULL)
    # dest_gfq = Base.unsafe_wrap(Vector{Int},pointer(pp.global_first_quadrant),MPI.Comm_size(MPI.COMM_WORLD)+1)
    return src_gfq, gfq, src_flt, src_llt
end
function get_receive_send(src_gfq::Vector, dest_gfq::Vector)
    is_receives = falses(MPI.Comm_size(MPI.COMM_WORLD))
    receive_nums = Vector{Int}(undef, MPI.Comm_size(MPI.COMM_WORLD))
    is_sends = falses(MPI.Comm_size(MPI.COMM_WORLD))
    send_nums = Vector{Int}(undef, MPI.Comm_size(MPI.COMM_WORLD))
    lnqs = [src_gfq[i+1] - src_gfq[i] for i = 1:MPI.Comm_size(MPI.COMM_WORLD)]
    nlnqs = [dest_gfq[i+1] - dest_gfq[i] for i = 1:MPI.Comm_size(MPI.COMM_WORLD)]
    lfq = dest_gfq[MPI.Comm_rank(MPI.COMM_WORLD)+1]
    llq = dest_gfq[MPI.Comm_rank(MPI.COMM_WORLD)+2] - 1
    r1 = MPI.Comm_rank(MPI.COMM_WORLD) == 0 ? 0 : findfirst(x -> x > lfq, src_gfq) - 2
    r2 =
        MPI.Comm_rank(MPI.COMM_WORLD) == MPI.Comm_size(MPI.COMM_WORLD) - 1 ?
        MPI.Comm_size(MPI.COMM_WORLD) - 1 : findfirst(x -> x > llq, src_gfq) - 2
    is_receives[r1+1:r2+1] .= true
    is_receives[MPI.Comm_rank(MPI.COMM_WORLD)+1] = false
    if r1 == r2
        receive_nums[r1+1] = llq - lfq + 1
    else
        for i in eachindex(is_receives)
            if is_receives[i]
                receive_nums[i] = min(src_gfq[i+1] - lfq, lnqs[i])
                if i == r2 + 1
                    receive_nums[i] = min(llq - src_gfq[i] + 1, llq - lfq + 1)
                end
            end
        end
    end
    lfq = src_gfq[MPI.Comm_rank(MPI.COMM_WORLD)+1]
    llq = src_gfq[MPI.Comm_rank(MPI.COMM_WORLD)+2] - 1
    s1 = MPI.Comm_rank(MPI.COMM_WORLD) == 0 ? 0 : findfirst(x -> x > lfq, dest_gfq) - 2
    s2 =
        MPI.Comm_rank(MPI.COMM_WORLD) == MPI.Comm_size(MPI.COMM_WORLD) - 1 ?
        MPI.Comm_size(MPI.COMM_WORLD) - 1 : findfirst(x -> x > llq, dest_gfq) - 2
    is_sends[s1+1:s2+1] .= true
    is_sends[MPI.Comm_rank(MPI.COMM_WORLD)+1] = false
    if s1 == s2
        send_nums[s1+1] = llq - lfq + 1
    else
        for i in eachindex(is_sends)
            if is_sends[i]
                send_nums[i] = min(dest_gfq[i+1] - lfq, nlnqs[i])
                if i == s2 + 1
                    send_nums[i] = min(llq - dest_gfq[i] + 1, llq - lfq + 1)
                end
            end
        end
    end
    receives = findall(is_receives)
    sends = findall(is_sends)
    receive_nums = receive_nums[is_receives]
    send_nums = send_nums[is_sends]
    return receives, sends, receive_nums, send_nums
end
function up_transfer_wrap(DIM::Integer,NDF::Integer,sends, send_nums, trees::Vector)
    s_datas = Vector{Transfer_Data{DIM,NDF}}(undef, 0)
    s_vs_numss = Vector{Vector{Int}}(undef, 0)
    send_index = 1
    index = 1
    up_num = length(sends)
    encs = Vector{Int}(undef, 0)
    ws = Vector{Cdouble}(undef, 0)
    vs_levels = Vector{Int8}(undef, 0)
    vs_midpoints = Vector{Cdouble}(undef, 0)
    vs_df = Vector{Cdouble}(undef, 0)
    s_vs_nums = Vector{Int}(undef, send_nums[send_index])
    for i in eachindex(trees)
        for j in eachindex(trees[i])
            sends[send_index] > MPI.Comm_rank(MPI.COMM_WORLD) + 1 && break
            ps_data = trees[i][j]
            if isa(ps_data,InsideSolidData)
                push!(encs,0)
                append!(encs,zeros(SOLID_CELL_ID_NUM))
                append!(ws,zeros(DIM+2))
                s_vs_nums[index] = 0
            else
                push!(encs,ps_data.bound_enc)
                append!(encs,ps_data.solid_cell_index)
                append!(ws, ps_data.w)
                append!(vs_levels, ps_data.vs_data.level)
                append!(vs_midpoints, reshape(ps_data.vs_data.midpoint, :))
                append!(vs_df, reshape(ps_data.vs_data.df, :))
                s_vs_nums[index] = ps_data.vs_data.vs_num
            end
            index += 1
            if index > send_nums[send_index]
                push!(s_datas, Transfer_Data{DIM,NDF}(encs, ws, vs_levels, vs_midpoints, vs_df))
                push!(s_vs_numss, s_vs_nums)
                index = 1
                send_index += 1
                send_index > up_num && break
                encs = Vector{Int}(undef, 0)
                ws = Vector{Cdouble}(undef, 0)
                vs_levels = Vector{Int8}(undef, 0)
                vs_midpoints = Vector{Cdouble}(undef, 0)
                vs_df = Vector{Cdouble}(undef, 0)
                s_vs_nums = Vector{Int}(undef, send_nums[send_index])
            end
        end
        send_index > up_num && break
        sends[send_index] > MPI.Comm_rank(MPI.COMM_WORLD) + 1 && break
    end
    return s_vs_numss, s_datas
end
function down_transfer_wrap(DIM::Integer,NDF::Integer,sends, send_nums, trees::Vector)
    s_datas = Vector{Transfer_Data{DIM,NDF}}(undef, 0)
    s_vs_numss = Vector{Vector{Int}}(undef, 0)
    send_index = length(sends)
    index = send_nums[send_index]
    encs = Vector{Int}(undef, 0)
    ws = Vector{Cdouble}(undef, 0)
    vs_levels = Vector{Int8}(undef, 0)
    vs_midpoints = Vector{Cdouble}(undef, 0)
    vs_df = Vector{Cdouble}(undef, 0)
    s_vs_nums = Vector{Int}(undef, send_nums[send_index])
    for i in reverse(eachindex(trees))
        for j in reverse(eachindex(trees[i]))
            sends[send_index] < MPI.Comm_rank(MPI.COMM_WORLD) + 1 && break
            ps_data = trees[i][j]
            if isa(ps_data,InsideSolidData)
                pushfirst!(encs,0)
                prepend!(encs,zeros(SOLID_CELL_ID_NUM))
                prepend!(ws,zeros(DIM+2))
                s_vs_nums[index] = 0
            else
                prepend!(encs,ps_data.solid_cell_index)
                pushfirst!(encs, ps_data.bound_enc)
                prepend!(ws, ps_data.w)
                prepend!(vs_levels, ps_data.vs_data.level)
                prepend!(vs_midpoints, reshape(ps_data.vs_data.midpoint, :))
                prepend!(vs_df, reshape(ps_data.vs_data.df, :))
                s_vs_nums[index] = ps_data.vs_data.vs_num
            end
            index -= 1
            if index < 1
                pushfirst!(s_datas, Transfer_Data{DIM,NDF}(encs, ws, vs_levels, vs_midpoints, vs_df))
                pushfirst!(s_vs_numss, s_vs_nums)
                send_index -= 1
                send_index < 1 && break
                index = send_nums[send_index]
                encs = Vector{Int}(undef, 0)
                ws = Vector{Cdouble}(undef, 0)
                vs_levels = Vector{Int8}(undef, 0)
                vs_midpoints = Vector{Cdouble}(undef, 0)
                vs_df = Vector{Cdouble}(undef, 0)
                s_vs_nums = Vector{Int}(undef, send_nums[send_index])
            end
        end
        send_index < 1 && break
        sends[send_index] < MPI.Comm_rank(MPI.COMM_WORLD) + 1 && break
    end
    return s_vs_numss, s_datas
end
function transfer_wrap(sends::Vector, send_nums::Vector, amr::AMR{DIM,NDF}) where{DIM,NDF}
    isempty(sends) && return Vector{Vector{Int}}(undef, 0), Vector{Transfer_Data{DIM,NDF}}(undef, 0)
    trees = amr.field.trees.data
    up_vs_numss, up_datas = up_transfer_wrap(DIM,NDF,sends, send_nums, trees)
    down_vs_numss, down_datas = down_transfer_wrap(DIM,NDF,sends, send_nums, trees)
    append!(up_vs_numss, down_vs_numss)
    append!(up_datas, down_datas)
    return up_vs_numss, up_datas
end
function up_transfer_pop!(sends, send_nums, amr::AMR)
    !any(x -> x < MPI.Comm_rank(MPI.COMM_WORLD) + 1, sends) && return nothing
    nums = 0
    trees = amr.field.trees.data
    for i in eachindex(sends)
        if sends[i] < MPI.Comm_rank(MPI.COMM_WORLD) + 1
            nums += send_nums[i]
        end
    end
    for i in eachindex(trees)
        for _ in eachindex(trees[i])
            popfirst!(trees[i])
            nums -= 1
            nums < 1 && break
        end
        nums < 1 && break
    end
end
function down_transfer_pop!(sends, send_nums, amr::AMR)
    !any(x -> x > MPI.Comm_rank(MPI.COMM_WORLD) + 1, sends) && return nothing
    nums = 0
    trees = amr.field.trees.data
    for i in eachindex(sends)
        if sends[i] > MPI.Comm_rank(MPI.COMM_WORLD) + 1
            nums += send_nums[i]
        end
    end
    for i in reverse(eachindex(trees))
        for _ in reverse(eachindex(trees[i]))
            pop!(trees[i])
            nums -= 1
            nums < 1 && break
        end
        nums < 1 && break
    end
end
function transfer_pop!(sends::Vector, send_nums::Vector, amr::AMR)
    up_transfer_pop!(sends, send_nums, amr)
    down_transfer_pop!(sends, send_nums, amr)
    trees = amr.field.trees.data
    deleteat!(trees, isempty.(trees))
end
function pre_transfer(
    sends::Vector{Int},
    s_vs_numss::Vector,
    receives::Vector{Int},
    receive_nums::Vector{Int},
)
    for i in eachindex(sends)
        sreq = MPI.Isend(
            s_vs_numss[i],
            MPI.COMM_WORLD;
            dest = sends[i] - 1,
            tag = COMM_NUMS_TAG + MPI.Comm_rank(MPI.COMM_WORLD),
        )
        MPI.Wait(sreq)
    end
    r_vs_numss = Vector{Vector{Int}}(undef, length(receives))
    for i in eachindex(receives)
        r_vs_nums = Vector{Int}(undef, receive_nums[i])
        rreq = MPI.Irecv!(
            r_vs_nums,
            MPI.COMM_WORLD;
            source = receives[i] - 1,
            tag = COMM_NUMS_TAG + receives[i] - 1,
        )
        MPI.Wait(rreq)
        r_vs_numss[i] = r_vs_nums
    end
    return r_vs_numss
end
function transfer(
    sends::Vector{Int},
    s_datas::Vector{Transfer_Data{DIM,NDF}},
    receives::Vector{Int},
    receive_nums::Vector{Int},
    r_vs_numss::Vector,
) where{DIM,NDF}
    reqs = Vector{MPI.Request}(undef, 0)
    for i in eachindex(sends)
        sreq = MPI.Isend(
            s_datas[i].encs,
            MPI.COMM_WORLD;
            dest = sends[i] - 1,
            tag = COMM_DATA_TAG + MPI.Comm_rank(MPI.COMM_WORLD),
        )
        push!(reqs, sreq)
        sreq = MPI.Isend(
            s_datas[i].w,
            MPI.COMM_WORLD;
            dest = sends[i] - 1,
            tag = COMM_DATA_TAG + MPI.Comm_rank(MPI.COMM_WORLD),
        )
        push!(reqs, sreq)
        sreq = MPI.Isend(
            s_datas[i].vs_levels,
            MPI.COMM_WORLD;
            dest = sends[i] - 1,
            tag = COMM_DATA_TAG + MPI.Comm_rank(MPI.COMM_WORLD),
        )
        push!(reqs, sreq)
        sreq = MPI.Isend(
            s_datas[i].vs_midpoints,
            MPI.COMM_WORLD;
            dest = sends[i] - 1,
            tag = COMM_DATA_TAG + MPI.Comm_rank(MPI.COMM_WORLD),
        )
        push!(reqs, sreq)
        sreq = MPI.Isend(
            s_datas[i].vs_df,
            MPI.COMM_WORLD;
            dest = sends[i] - 1,
            tag = COMM_DATA_TAG + MPI.Comm_rank(MPI.COMM_WORLD),
        )
        push!(reqs, sreq)
    end
    r_datas = Vector{Transfer_Data{DIM,NDF}}(undef, 0)
    for i in eachindex(receives)
        total_vs_num = sum(r_vs_numss[i])
        r_data = Transfer_Data(DIM,NDF,receive_nums[i], total_vs_num)
        rreq = MPI.Irecv!(
            r_data.encs,
            MPI.COMM_WORLD;
            source = receives[i] - 1,
            tag = COMM_DATA_TAG + receives[i] - 1,
        )
        push!(reqs, rreq)
        rreq = MPI.Irecv!(
            r_data.w,
            MPI.COMM_WORLD;
            source = receives[i] - 1,
            tag = COMM_DATA_TAG + receives[i] - 1,
        )
        push!(reqs, rreq)
        rreq = MPI.Irecv!(
            r_data.vs_levels,
            MPI.COMM_WORLD;
            source = receives[i] - 1,
            tag = COMM_DATA_TAG + receives[i] - 1,
        )
        push!(reqs, rreq)
        rreq = MPI.Irecv!(
            r_data.vs_midpoints,
            MPI.COMM_WORLD;
            source = receives[i] - 1,
            tag = COMM_DATA_TAG + receives[i] - 1,
        )
        push!(reqs, rreq)
        rreq = MPI.Irecv!(
            r_data.vs_df,
            MPI.COMM_WORLD;
            source = receives[i] - 1,
            tag = COMM_DATA_TAG + receives[i] - 1,
        )
        push!(reqs, rreq)
        push!(r_datas, r_data)
    end
    MPI.Waitall!(reqs)
    up_index = [x < MPI.Comm_rank(MPI.COMM_WORLD) + 1 for x in receives]
    down_index = [x > MPI.Comm_rank(MPI.COMM_WORLD) + 1 for x in receives]
    up_vs_numss = r_vs_numss[up_index]
    down_vs_numss = r_vs_numss[down_index]
    up_datas = r_datas[up_index]
    down_datas = r_datas[down_index]
    return up_vs_numss, down_vs_numss, up_datas, down_datas
end
function partition_transfer(
    sends::Vector{Int},
    s_vs_numss::Vector,
    s_datas::Vector{Transfer_Data{DIM,NDF}},
    receives::Vector{Int},
    receive_nums::Vector{Int},
)where{DIM,NDF}
    r_vs_numss = pre_transfer(sends, s_vs_numss, receives, receive_nums)
    transfer(sends, s_datas, receives, receive_nums, r_vs_numss)
end
function unpack_data(vs_nums, data, amr::AMR{DIM,NDF}) where{DIM,NDF}
    transfer_ps_datas = Array{AbstractPsData{DIM,NDF}}(undef, length(vs_nums))
    quadrature = amr.global_data.config.quadrature
    vs_trees_num = reduce(*, amr.global_data.config.vs_trees_num)
    vs_space = 1.0
    for i = 1:DIM
        vs_space *= quadrature[2*i] - quadrature[2*i-1]
    end
    tree_weight = vs_space / vs_trees_num
    offset = 0
    for i in eachindex(vs_nums)
        vs_num = vs_nums[i]
        if vs_num==0
            transfer_ps_datas[i] = InsideSolidData{DIM,NDF}()
        else
            encs = data.encs[(SOLID_CELL_ID_NUM+1)*(i-1)+1:(SOLID_CELL_ID_NUM+1)*i]
            w = data.w[(DIM+2)*(i-1)+1:(DIM+2)*i]
            vs_levels = data.vs_levels[offset+1:offset+vs_num]
            vs_midpoints =
                reshape(data.vs_midpoints[DIM*offset+1:DIM*(offset+vs_num)], vs_num, DIM)
            vs_df = reshape(data.vs_df[NDF*offset+1:NDF*(offset+vs_num)], vs_num, NDF)
            vs_weight = @. tree_weight / 2.0^(DIM * vs_levels)
            vs_data = VS_Data{DIM,NDF}(
                vs_num,
                vs_levels,
                vs_weight,
                vs_midpoints,
                vs_df,
                zeros(vs_num, NDF, DIM),
                zeros(vs_num, NDF),
            )
            transfer_ps_datas[i] = PS_Data(DIM,NDF,encs,w, vs_data)
        end
        offset += vs_num
    end
    return transfer_ps_datas
end
function unpack_data(up_vs_numss, down_vs_numss, up_datas, down_datas, amr::AMR{DIM,NDF}) where{DIM,NDF}
    up_ps_datas = Vector{AbstractPsData{DIM,NDF}}(undef, 0)
    down_ps_datas = Vector{AbstractPsData{DIM,NDF}}(undef, 0)
    for i in eachindex(up_vs_numss)
        append!(up_ps_datas, unpack_data(up_vs_numss[i], up_datas[i], amr))
    end
    for i in eachindex(down_vs_numss)
        append!(down_ps_datas, unpack_data(down_vs_numss[i], down_datas[i], amr))
    end
    return up_ps_datas, down_ps_datas
end
function init_up_quadrants!(ip, dp, ti_data::Transfer_Init, treeid::Integer, tree_datas)
    ti_data.up_index > ti_data.up_num && return nothing
    ps_data = ti_data.up_data[ti_data.up_index]
    if !isa(ps_data,InsideSolidData)
        ps_data.ds, ps_data.midpoint = quad_to_cell(ip.p4est, treeid, ip.quad)
    end
    dp[] = P4est_PS_Data(pointer_from_objref(ps_data))
    if treeid < ti_data.old_flt
        push!(tree_datas, ps_data)
    else
        insert!(tree_datas, ti_data.up_insert_index, ps_data)
        ti_data.up_insert_index += 1
    end
    ti_data.up_index += 1
end
function init_down_quadrants!(ip, dp, ti_data::Transfer_Init, treeid::Integer, tree_datas)
    local_quadid(ip) < (ip.p4est.local_num_quadrants[] - ti_data.down_num) &&
        return nothing
    ps_data = ti_data.down_data[ti_data.down_index]
    if !isa(ps_data,InsideSolidData)
        ps_data.ds, ps_data.midpoint = quad_to_cell(ip.p4est, treeid, ip.quad)
    end
    dp[] = P4est_PS_Data(pointer_from_objref(ps_data))
    push!(tree_datas, ps_data)
    ti_data.down_index += 1
end
function init_transferred_quadrant!(ip, data, dp)
    ti_data = unsafe_pointer_to_objref(data)
    trees = ti_data.amr.field.trees
    treeid = ip.treeid[] - trees.offset
    tree_datas = trees.data[treeid]
    init_up_quadrants!(ip, dp, ti_data, ip.treeid[], tree_datas)
    init_down_quadrants!(ip, dp, ti_data, ip.treeid[], tree_datas)
end
function init_transferred_quadrant!(info, data)
    AMR_volume_iterate(info, data, P4est_PS_Data, init_transferred_quadrant!)
end
function init_transferred_quadrant!(p4est::P_pxest_t, ti_data::Transfer_Init)
    ghost = ti_data.amr.global_data.forest.ghost
    p_data = pointer_from_objref(ti_data)
    GC.@preserve ti_data AMR_4est_volume_iterate(
        p4est,
        ghost,
        p_data,
        init_transferred_quadrant!,
    )
end
function init_transferred_ps!(p4est::P_pxest_t, ti_data::Transfer_Init)
    insert_trees!(p4est, ti_data.amr, ti_data)
    init_transferred_quadrant!(p4est, ti_data)
end
function insert_trees!(p4est::P_pxest_t, amr::AMR{DIM,NDF}, ti_data::Transfer_Init) where{DIM,NDF}
    trees = amr.field.trees.data
    pp = PointerWrapper(p4est)
    if !isempty(trees)
        if pp.first_local_tree[] < ti_data.old_flt
            for i = 1:ti_data.old_flt-pp.first_local_tree[]
                pushfirst!(trees, Vector{PS_Data{DIM,NDF}}(undef, 0))
            end
        end
        if pp.last_local_tree[] > ti_data.old_llt
            for i = 1:pp.last_local_tree[]-ti_data.old_llt
                push!(trees, Vector{PS_Data{DIM,NDF}}(undef, 0))
            end
        end
    else
        for _ = 1:pp.last_local_tree[]-pp.first_local_tree[]+1
            push!(trees, Vector{PS_Data{DIM,NDF}}(undef, 0))
        end
    end
    amr.field.trees.offset = pp.first_local_tree[] - 1
end

function ps_partition!(p4est::P_pxest_t, amr::AMR)
    src_gfq, dest_gfq, src_flt, src_llt = partition!(p4est)
    receives, sends, receive_nums, send_nums = get_receive_send(src_gfq, dest_gfq)
    s_vs_numss, s_datas = transfer_wrap(sends, send_nums, amr)
    transfer_pop!(sends, send_nums, amr)
    up_vs_numss, down_vs_numss, up_datas, down_datas =
        partition_transfer(sends, s_vs_numss, s_datas, receives, receive_nums)
    up_ps_datas, down_ps_datas =
        unpack_data(up_vs_numss, down_vs_numss, up_datas, down_datas, amr)
    ti_data = Transfer_Init(
        length(up_ps_datas),
        length(down_ps_datas),
        up_ps_datas,
        down_ps_datas,
        src_flt,
        src_llt,
        1,
        1,
        1,
        amr,
    )
    init_transferred_ps!(p4est, ti_data)
end
