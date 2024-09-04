size_Ghost_Data(vs_num,DIM,NDF) = 3 * DIM + 3 + NDF * vs_num
size_Ghost_Slope(vs_num,DIM,NDF) = vs_num * DIM * NDF
size_Ghost_VS_Structure(vs_num,DIM) = vs_num * (DIM + 2)
function get_vs_num(forest::P_pxest_t, ghost::P_pxest_ghost_t)
    get_vs_num(PointerWrapper(forest), PointerWrapper(ghost))
end
function pw_mirror_quadrant(pp::PointerWrapper{p4est_t},gp::PointerWrapper{p4est_ghost_t},i::Integer)
        pm = iPointerWrapper(gp.mirrors, p4est_quadrant_t, i - 1)
        pt = iPointerWrapper(pp.trees, p4est_tree_t, pm.p.piggy3.which_tree[])
        pq = iPointerWrapper(
            pt.quadrants,
            p4est_quadrant_t,
            pm.p.piggy3.local_num[] - pt.quadrants_offset[],
        )
end
function pw_mirror_quadrant(pp::PointerWrapper{p8est_t},gp::PointerWrapper{p8est_ghost_t},i::Integer)
    pm = iPointerWrapper(gp.mirrors, p8est_quadrant_t, i - 1)
    pt = iPointerWrapper(pp.trees, p8est_tree_t, pm.p.piggy3.which_tree[])
    pq = iPointerWrapper(
        pt.quadrants,
        p8est_quadrant_t,
        pm.p.piggy3.local_num[] - pt.quadrants_offset[],
    )
end
function get_vs_num(pp::PW_pxest_t, gp::PW_pxest_ghost_t)
    vs_num = 0
    for i = 1:gp.mirrors.elem_count[]
        pq = pw_mirror_quadrant(pp,gp,i)
        dp = PointerWrapper(P4est_PS_Data, pq.p.user_data[])
        ps_data = unsafe_pointer_to_objref(pointer(dp.ps_data))
        vs_num = max(vs_num, ps_data.vs_data.vs_num)
    end
    vs_num = MPI.Allreduce(vs_num, max, MPI.COMM_WORLD)
end
function get_mirror_data_inner!(ps_data::PS_Data{DIM,NDF}, vs_temp::AbstractVector) where{DIM,NDF}
    vs_data = ps_data.vs_data
    vs_num = vs_data.vs_num
    vs_temp[1:NDF*vs_num] .= reshape(vs_data.df, vs_num * NDF)
end
function get_mirror_data(ps4est, global_data::Global_Data{DIM,NDF}) where{DIM,NDF}
    GC.@preserve ps4est global_data begin
        gp = PointerWrapper(global_data.forest.ghost)
        pp = PointerWrapper(ps4est)
        global_data.status.max_vs_num = vs_num = get_vs_num(pp, gp)
        mirror_data_pointers = Array{Ptr{Cdouble}}(undef, gp.mirrors.elem_count[])
        for i = 1:gp.mirrors.elem_count[]
            pq = pw_mirror_quadrant(pp,gp,i)
            dp = PointerWrapper(P4est_PS_Data, pq.p.user_data[])
            ps_data = unsafe_pointer_to_objref(pointer(dp.ps_data))
            p = Ptr{Cdouble}(sc_malloc(-1, (3 * DIM + 3 + NDF * vs_num) * sizeof(Cdouble)))
            ap = Base.unsafe_wrap(Vector{Cdouble}, p, 3 * DIM + 3 + NDF * vs_num)
            offset = 0
            ap[1:(offset+=DIM)] .= ps_data.ds
            ap[offset+1:(offset+=DIM)] .= ps_data.midpoint
            ap[offset+1:(offset+=DIM+2)] .= ps_data.w
            ap[offset+=1] = ps_data.vs_data.vs_num
            vs_temp = @view(ap[offset+1:(offset+=NDF*vs_num)])
            get_mirror_data_inner!(ps_data, vs_temp)
            mirror_data_pointers[i] = p
        end
    end
    return mirror_data_pointers
end
function get_mirror_slope_inner!(ps_data::PS_Data{DIM,NDF}, vs_temp::AbstractVector) where{DIM,NDF}
    vs_data = ps_data.vs_data
    vs_num = size(vs_data.weight, 1)
    vs_temp[1:(vs_num*NDF*DIM)] .= reshape(vs_data.sdf, vs_num * NDF * DIM)
end
function get_mirror_slope(ps4est, global_data::Global_Data{DIM,NDF}) where{DIM,NDF}
    GC.@preserve ps4est global_data begin
        gp = PointerWrapper(global_data.forest.ghost)
        pp = PointerWrapper(ps4est)
        vs_num = global_data.status.max_vs_num
        mirror_slope_pointers = Array{Ptr{Cdouble}}(undef, gp.mirrors.elem_count[])
        for i = 1:gp.mirrors.elem_count[]
            pq = pw_mirror_quadrant(pp,gp,i)
            dp = PointerWrapper(P4est_PS_Data, pq.p.user_data[])
            ps_data = unsafe_pointer_to_objref(pointer(dp.ps_data))
            p = Ptr{Cdouble}(sc_malloc(-1, (vs_num * DIM * NDF) * sizeof(Cdouble)))
            ap = Base.unsafe_wrap(Vector{Cdouble}, p, vs_num * DIM * NDF)
            get_mirror_slope_inner!(ps_data, ap)
            mirror_slope_pointers[i] = p
        end
    end
    return mirror_slope_pointers
end
function get_mirror_structure_inner!(ps_data::PS_Data{DIM}, weight_temp, level_temp, midpoint_temp) where{DIM}
    vs_data = ps_data.vs_data
    vs_num = vs_data.vs_num
    weight_temp[1:vs_num] .= vs_data.weight
    level_temp[1:vs_num] .= vs_data.level
    midpoint_temp[1:(vs_num*DIM)] .= reshape(vs_data.midpoint, vs_num * DIM)
end
function get_mirror_structure(ps4est, global_data::Global_Data{DIM,NDF}) where{DIM,NDF}
    GC.@preserve ps4est global_data begin
        gp = PointerWrapper(global_data.forest.ghost)
        pp = PointerWrapper(ps4est)
        vs_num = global_data.status.max_vs_num
        mirror_structure_pointers = Array{Ptr{Cdouble}}(undef, gp.mirrors.elem_count[])
        for i = 1:gp.mirrors.elem_count[]
            pq = pw_mirror_quadrant(pp,gp,i)
            dp = PointerWrapper(P4est_PS_Data, pq.p.user_data[])
            ps_data = unsafe_pointer_to_objref(pointer(dp.ps_data))
            p = Ptr{Cdouble}(sc_malloc(-1, vs_num * (DIM + 2) * sizeof(Cdouble)))
            ap = Base.unsafe_wrap(Vector{Cdouble}, p, vs_num * (DIM + 2))
            weight_temp = @view(ap[1:vs_num])
            level_temp = @view(ap[vs_num+1:vs_num*2])
            midpoint_temp = @view(ap[vs_num*2+1:vs_num*(DIM+2)])
            get_mirror_structure_inner!(ps_data, weight_temp, level_temp, midpoint_temp)
            mirror_structure_pointers[i] = p
        end
        return mirror_structure_pointers
    end
end
function update_mirror_data!(ps4est, amr::AMR{DIM,NDF}) where{DIM,NDF}
    GC.@preserve ps4est amr begin
        mirror_data_pointers = amr.ghost.ghost_exchange.mirror_data_pointers
        vs_num = amr.global_data.status.max_vs_num
        gp = PointerWrapper(amr.global_data.forest.ghost)
        pp = PointerWrapper(ps4est)
        for i in eachindex(mirror_data_pointers)
            pq = pw_mirror_quadrant(pp,gp,i)
            dp = PointerWrapper(P4est_PS_Data, pq.p.user_data[])
            ps_data = unsafe_pointer_to_objref(pointer(dp.ps_data))
            ap = Base.unsafe_wrap(
                Vector{Cdouble},
                mirror_data_pointers[i],
                3 * DIM + 3 + NDF * vs_num,
            )
            vs_temp = @view(ap[3*DIM+3+1:end])
            get_mirror_data_inner!(ps_data, vs_temp)
        end
    end
end
function update_mirror_slope!(ps4est, amr::AMR{DIM,NDF}) where{DIM,NDF}
    GC.@preserve ps4est amr begin
        mirror_slope_pointers = amr.ghost.ghost_exchange.mirror_slope_pointers
        vs_num = amr.global_data.status.max_vs_num
        gp = PointerWrapper(amr.global_data.forest.ghost)
        pp = PointerWrapper(ps4est)
        for i in eachindex(mirror_slope_pointers)
            pq = pw_mirror_quadrant(pp,gp,i)
            dp = PointerWrapper(P4est_PS_Data, pq.p.user_data[])
            ps_data = unsafe_pointer_to_objref(pointer(dp.ps_data))
            ap = Base.unsafe_wrap(
                Vector{Cdouble},
                mirror_slope_pointers[i],
                vs_num * DIM * NDF,
            )
            get_mirror_slope_inner!(ps_data, ap)
        end
    end
end
function ghost_data_alloc(N::Int, ghost::P_pxest_ghost_t)
    Ptr{Cdouble}(sc_malloc(-1, sizeof(Cdouble) * N * PointerWrapper(ghost).ghosts.elem_count[]))
end
function ghost_data_alloc(::Type{T},N::Int, ghost::P_pxest_ghost_t) where{T}
    Ptr{T}(sc_malloc(-1, sizeof(Cdouble) * N * PointerWrapper(ghost).ghosts.elem_count[]))
end
function amr_exchange(
    p4est::Ptr{p4est_t},
    ghost::Ptr{p4est_ghost_t},
    mirror_data_pointers::Array{Ptr{Cdouble}},
    N::Int,
)
    GC.@preserve mirror_data_pointers ghost N begin
        ghost_datas = ghost_data_alloc(N, ghost)
        GC.@preserve ghost_datas p4est_ghost_exchange_custom(
            p4est,
            ghost,
            sizeof(Cdouble) * N,
            mirror_data_pointers,
            ghost_datas,
        )
    end
    return ghost_datas
end
function amr_exchange!(
    p4est::Ptr{p4est_t},
    ghost::Ptr{p4est_ghost_t},
    ghost_datas::Ptr{Cdouble},
    mirror_data_pointers::Array{Ptr{Cdouble}},
    N::Int,
)
    GC.@preserve mirror_data_pointers ghost ghost_datas begin
        p4est_ghost_exchange_custom(
            p4est,
            ghost,
            sizeof(Cdouble) * N,
            mirror_data_pointers,
            ghost_datas,
        )
    end
end
function amr_exchange(
    p4est::Ptr{p8est_t},
    ghost::Ptr{p8est_ghost_t},
    mirror_data_pointers::Array{Ptr{Cdouble}},
    N::Int,
)
    GC.@preserve mirror_data_pointers ghost N begin
        ghost_datas = ghost_data_alloc(N, ghost)
        GC.@preserve ghost_datas p8est_ghost_exchange_custom(
            p4est,
            ghost,
            sizeof(Cdouble) * N,
            mirror_data_pointers,
            ghost_datas,
        )
    end
    return ghost_datas
end
function amr_exchange!(
    p4est::Ptr{p8est_t},
    ghost::Ptr{p8est_ghost_t},
    ghost_datas::Ptr{Cdouble},
    mirror_data_pointers::Array{Ptr{Cdouble}},
    N::Int,
)
    GC.@preserve mirror_data_pointers ghost ghost_datas begin
        p8est_ghost_exchange_custom(
            p4est,
            ghost,
            sizeof(Cdouble) * N,
            mirror_data_pointers,
            ghost_datas,
        )
    end
end
function data_exchange!(p4est::P_pxest_t, amr::AMR{DIM,NDF}) where{DIM,NDF}
    MPI.Comm_size(MPI.COMM_WORLD) == 1 && return nothing
    update_mirror_data!(p4est, amr)
    amr_exchange!(
        p4est,
        amr.global_data.forest.ghost,
        amr.ghost.ghost_exchange.ghost_datas,
        amr.ghost.ghost_exchange.mirror_data_pointers,
        size_Ghost_Data(amr.global_data.status.max_vs_num,DIM,NDF),
    )
end
function slope_exchange!(p4est::P_pxest_t, amr::AMR{DIM,NDF}) where{DIM,NDF}
    MPI.Comm_size(MPI.COMM_WORLD) == 1 && return nothing
    update_mirror_slope!(p4est, amr)
    amr_exchange!(
        p4est,
        amr.global_data.forest.ghost,
        amr.ghost.ghost_exchange.ghost_slopes,
        amr.ghost.ghost_exchange.mirror_slope_pointers,
        size_Ghost_Slope(amr.global_data.status.max_vs_num,DIM,NDF),
    )
end
function initialize_ghost_exchange(ps4est, global_data::Global_Data{DIM,NDF}) where{DIM,NDF}
    mirror_data_pointers = get_mirror_data(ps4est, global_data)
    ghost_datas = amr_exchange(
        ps4est,
        global_data.forest.ghost,
        mirror_data_pointers,
        size_Ghost_Data(global_data.status.max_vs_num,DIM,NDF),
    )
    mirror_slope_pointers = get_mirror_slope(ps4est, global_data)
    ghost_slopes = amr_exchange(
        ps4est,
        global_data.forest.ghost,
        mirror_slope_pointers,
        size_Ghost_Slope(global_data.status.max_vs_num,DIM,NDF),
    )
    mirror_structure_pointers = get_mirror_structure(ps4est, global_data)
    ghost_structures = amr_exchange(
        ps4est,
        global_data.forest.ghost,
        mirror_structure_pointers,
        size_Ghost_VS_Structure(global_data.status.max_vs_num,DIM),
    )
    return Ghost_Exchange(
        ghost_datas,
        ghost_slopes,
        ghost_structures,
        mirror_data_pointers,
        mirror_slope_pointers,
        mirror_structure_pointers,
    )
end

function ghost_data_ptr(N::Int, ghost_data::Ptr{T}, qid::Integer) where {T}
    return Ptr{Cdouble}(ghost_data + qid * sizeof(T) * N)
end
function initialize_ghost_wrap(global_data::Global_Data{DIM,NDF}, ghost_exchange::Ghost_Exchange) where{DIM,NDF}
    ghost_wrap =
        Array{Ghost_PS_Data}(undef, PointerWrapper(global_data.forest.ghost).ghosts.elem_count[])
    global_vs_num = global_data.status.max_vs_num
    for i in eachindex(ghost_wrap)
        p = ghost_data_ptr(
            size_Ghost_Data(global_vs_num,DIM,NDF),
            ghost_exchange.ghost_datas,
            i - 1,
        )
        offset = 0
        ds = Base.unsafe_wrap(Vector{Cdouble}, p + offset * sizeof(Cdouble), DIM)
        offset += DIM
        midpoint = Base.unsafe_wrap(Vector{Cdouble}, p + offset * sizeof(Cdouble), DIM)
        offset += DIM
        w = Base.unsafe_wrap(Vector{Cdouble}, p + offset * sizeof(Cdouble), DIM + 2)
        offset += DIM + 2
        vs_num = Int(Base.unsafe_load(p + offset * sizeof(Cdouble)))
        offset += 1
        # midpoint_vs = Base.unsafe_wrap(Matrix{Cdouble},Ptr{Cdouble}(pointer(ghost_data.micro)),(vs_num,DIM))
        df = Base.unsafe_wrap(Matrix{Cdouble}, p + offset * sizeof(Cdouble), (vs_num, NDF))
        p = ghost_data_ptr(
            size_Ghost_Slope(global_vs_num,DIM,NDF),
            ghost_exchange.ghost_slopes,
            i - 1,
        )
        sdf = Base.unsafe_wrap(Array{Cdouble}, p, (vs_num, NDF, DIM))
        p = ghost_data_ptr(
            size_Ghost_VS_Structure(global_vs_num,DIM),
            ghost_exchange.ghost_structures,
            i - 1,
        )
        offset = 0
        weight = Base.unsafe_wrap(Vector{Cdouble}, p, vs_num)
        offset += global_vs_num
        level = Base.unsafe_wrap(Vector{Cdouble}, p + offset * sizeof(Cdouble), vs_num)
        offset += global_vs_num
        midpoint_vs =
            Base.unsafe_wrap(Matrix{Cdouble}, p + offset * sizeof(Cdouble), (vs_num, DIM))
        vs_data = Ghost_VS_Data{DIM,NDF}(vs_num, Int.(level), weight, midpoint_vs, df, sdf)
        ghost_wrap[i] = Ghost_PS_Data{DIM,NDF}(ds, midpoint, w, vs_data)
    end
    return ghost_wrap
end

function update_ghost!(p4est::Ptr{p4est_t}, amr::AMR)
    ghost_exchange = amr.ghost.ghost_exchange
    global_data = amr.global_data
    finalize_ghost!(ghost_exchange)
    p4est_ghost_destroy(global_data.forest.ghost)
    global_data.forest.ghost = p8est_ghost_new(p4est, P4EST_CONNECT_FACE)
    amr.ghost.ghost_exchange = initialize_ghost_exchange(p4est, global_data)
    amr.ghost.ghost_wrap = initialize_ghost_wrap(global_data, amr.ghost.ghost_exchange)
end
function update_ghost!(p4est::Ptr{p8est_t}, amr::AMR)
    ghost_exchange = amr.ghost.ghost_exchange
    global_data = amr.global_data
    finalize_ghost!(ghost_exchange)
    p8est_ghost_destroy(global_data.forest.ghost)
    global_data.forest.ghost = p8est_ghost_new(p4est, P8EST_CONNECT_FACE)
    amr.ghost.ghost_exchange = initialize_ghost_exchange(p4est, global_data)
    amr.ghost.ghost_wrap = initialize_ghost_wrap(global_data, amr.ghost.ghost_exchange)
end
