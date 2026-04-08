size_Ghost_Data(vs_num,DIM,NDF) = 3 * DIM + 4 + NDF * vs_num
size_Ghost_Slope(vs_num,DIM,NDF) = vs_num * DIM * NDF+DIM*(DIM+2)
size_Ghost_VS_Structure(vs_num,DIM) = vs_num * (DIM + 2)
"""
$(TYPEDSIGNATURES)
Get the globally largest number of the velocity cells in ghost layers.
"""
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
        dp = PointerWrapper(P4estPsData, pq.p.user_data[])
        ps_data = unsafe_pointer_to_objref(pointer(dp.ps_data))
        isa(ps_data,InsideSolidData)&&continue
        vs_num = max(vs_num, ps_data.vs_data.vs_num)
    end
    vs_num = MPI.Allreduce(vs_num, max, MPI.COMM_WORLD)
end
function get_mirror_data_inner!(ps_data::PsData{DIM,NDF}, vs_temp::AbstractVector) where{DIM,NDF}
    vs_data = ps_data.vs_data
    vs_num = vs_data.vs_num
    vs_temp[1:NDF*vs_num] .= reshape(vs_data.df, vs_num * NDF)
end
# function get_mirror_data_inner!(::InsideSolidData,vs_temp::AbstractVector)
#     vs_temp[1] = EPS # Flag indicating InsideSolidData
# end
function get_mirror_data(p4est, kinfo::KInfo{DIM,NDF}) where{DIM,NDF}
    GC.@preserve p4est kinfo begin
        gp = PointerWrapper(kinfo.forest.ghost)
        pp = PointerWrapper(p4est)
        kinfo.status.max_vs_num = vs_num = get_vs_num(pp, gp)
        mirror_data_pointers = Array{Ptr{Cdouble}}(undef, gp.mirrors.elem_count[])
        for i = 1:gp.mirrors.elem_count[]
            pq = pw_mirror_quadrant(pp,gp,i)
            dp = PointerWrapper(P4estPsData, pq.p.user_data[])
            ps_data = unsafe_pointer_to_objref(pointer(dp.ps_data))
            p = Ptr{Cdouble}(sc_malloc(P4est.package_id(), (3 * DIM + 4 + NDF * vs_num) * sizeof(Cdouble)))
            ap = Base.unsafe_wrap(Vector{Cdouble}, p, 3 * DIM + 4 + NDF * vs_num)
            offset = 0
            if isa(ps_data,InsideSolidData)
                ap[1:(offset+=DIM)] .= ps_data.ds
                ap[offset+1:(offset+=DIM)] .= ps_data.midpoint
                ap[offset+1:(offset+=DIM+2)] .= Inf
                ap[offset+=1] = 0
                ap[offset+=1] = ps_data.bound_enc
            else
                ap[1:(offset+=DIM)] .= ps_data.ds
                ap[offset+1:(offset+=DIM)] .= ps_data.midpoint
                ap[offset+1:(offset+=DIM+2)] .= ps_data.w
                ap[offset+=1] = ps_data.vs_data.vs_num
                ap[offset+=1] = ps_data.bound_enc
                vs_temp = @view(ap[offset+1:(offset+=NDF*vs_num)])
                get_mirror_data_inner!(ps_data, vs_temp)
            end
            mirror_data_pointers[i] = p
        end
    end
    return mirror_data_pointers
end
function get_mirror_slope_inner!(ps_data::PsData{DIM,NDF}, vs_temp::AbstractVector) where{DIM,NDF}
    vs_data = ps_data.vs_data
    vs_num = size(vs_data.weight, 1)
    vs_temp[1:(vs_num*NDF*DIM)] .= reshape(vs_data.sdf, vs_num * NDF * DIM)
    vs_temp[vs_num*NDF*DIM+1:vs_num*NDF*DIM+DIM*(DIM+2)] .= reshape(ps_data.sw, DIM*(DIM+2))
end
function get_mirror_slope_inner!(::InsideSolidData, vs_temp::AbstractVector) 
    vs_temp[1] = EPS
end
function get_mirror_slope(p4est, kinfo::KInfo{DIM,NDF}) where{DIM,NDF}
    GC.@preserve p4est kinfo begin
        gp = PointerWrapper(kinfo.forest.ghost)
        pp = PointerWrapper(p4est)
        vs_num = kinfo.status.max_vs_num
        mirror_slope_pointers = Array{Ptr{Cdouble}}(undef, gp.mirrors.elem_count[])
        for i = 1:gp.mirrors.elem_count[]
            pq = pw_mirror_quadrant(pp,gp,i)
            dp = PointerWrapper(P4estPsData, pq.p.user_data[])
            ps_data = unsafe_pointer_to_objref(pointer(dp.ps_data))
            p = Ptr{Cdouble}(sc_malloc(P4est.package_id(), (vs_num * DIM * NDF+DIM*(DIM+2)) * sizeof(Cdouble)))
            ap = Base.unsafe_wrap(Vector{Cdouble}, p, vs_num * DIM * NDF+DIM*(DIM+2))
            get_mirror_slope_inner!(ps_data, ap)
            mirror_slope_pointers[i] = p
        end
    end
    return mirror_slope_pointers
end
function get_mirror_structure_inner!(ps_data::PsData{DIM}, weight_temp, level_temp, midpoint_temp) where{DIM}
    vs_data = ps_data.vs_data
    vs_num = vs_data.vs_num
    weight_temp[1:vs_num] .= vs_data.weight
    level_temp[1:vs_num] .= vs_data.level
    midpoint_temp[1:(vs_num*DIM)] .= reshape(vs_data.midpoint, vs_num * DIM)
end
function get_mirror_structure_inner!(ps_data::InsideSolidData, weight_temp, level_temp, midpoint_temp)
    weight_temp[1] = EPS
    level_temp[1] = EPS
    midpoint_temp[1] = EPS
end
function get_mirror_structure(p4est, kinfo::KInfo{DIM,NDF}) where{DIM,NDF}
    GC.@preserve p4est kinfo begin
        gp = PointerWrapper(kinfo.forest.ghost)
        pp = PointerWrapper(p4est)
        vs_num = kinfo.status.max_vs_num
        mirror_structure_pointers = Array{Ptr{Cdouble}}(undef, gp.mirrors.elem_count[])
        for i = 1:gp.mirrors.elem_count[]
            pq = pw_mirror_quadrant(pp,gp,i)
            dp = PointerWrapper(P4estPsData, pq.p.user_data[])
            ps_data = unsafe_pointer_to_objref(pointer(dp.ps_data))
            p = Ptr{Cdouble}(sc_malloc(P4est.package_id(), vs_num * (DIM + 2) * sizeof(Cdouble)))
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
function update_mirror_data!(p4est, ka::KA{DIM,NDF}) where{DIM,NDF}
    GC.@preserve p4est ka begin
        mirror_data_pointers = ka.kdata.ghost.ghost_pointers.mirror_data_pointers
        vs_num = ka.kinfo.status.max_vs_num
        gp = PointerWrapper(ka.kinfo.forest.ghost)
        pp = PointerWrapper(p4est)
        for i in eachindex(mirror_data_pointers)
            pq = pw_mirror_quadrant(pp,gp,i)
            dp = PointerWrapper(P4estPsData, pq.p.user_data[])
            ps_data = unsafe_pointer_to_objref(pointer(dp.ps_data))
            isa(ps_data,InsideSolidData)&&continue
            ap = Base.unsafe_wrap(
                Vector{Cdouble},
                mirror_data_pointers[i],
                3 * DIM + 4 + NDF * vs_num,
            )
			ap[DIM*2+1:DIM*3+2] .= ps_data.w
            vs_temp = @view(ap[3*DIM+4+1:end])
            get_mirror_data_inner!(ps_data, vs_temp)
        end
    end
end
function update_solid_mirror_data!(p4est,ka::KA{DIM,NDF}) where{DIM,NDF}
    GC.@preserve p4est ka begin
        mirror_data_pointers = ka.kdata.ghost.ghost_pointers.mirror_data_pointers
        vs_num = ka.kinfo.status.max_vs_num
        gp = PointerWrapper(ka.kinfo.forest.ghost)
        pp = PointerWrapper(p4est)
        for i in eachindex(mirror_data_pointers)
            pq = pw_mirror_quadrant(pp,gp,i)
            dp = PointerWrapper(P4estPsData, pq.p.user_data[])
            ps_data = unsafe_pointer_to_objref(pointer(dp.ps_data))
            (isa(ps_data,InsideSolidData)||ps_data.bound_enc>=0)&&continue
            ap = Base.unsafe_wrap(
                Vector{Cdouble},
                mirror_data_pointers[i],
                3 * DIM + 4 + NDF * vs_num, # ds(DIM), midpoint(DIM), w(DIM+2), vs_num(1), bound_enc(1)
            )
			ap[DIM*2+1:DIM*3+2] .= ps_data.w
            vs_temp = @view(ap[3*DIM+4+1:end])
            get_mirror_data_inner!(ps_data, vs_temp)
        end
    end
end
function update_mirror_slope!(p4est, ka::KA{DIM,NDF}) where{DIM,NDF}
    GC.@preserve p4est ka begin
        mirror_slope_pointers = ka.kdata.ghost.ghost_pointers.mirror_slope_pointers
        vs_num = ka.kinfo.status.max_vs_num
        gp = PointerWrapper(ka.kinfo.forest.ghost)
        pp = PointerWrapper(p4est)
        for i in eachindex(mirror_slope_pointers)
            pq = pw_mirror_quadrant(pp,gp,i)
            dp = PointerWrapper(P4estPsData, pq.p.user_data[])
            ps_data = unsafe_pointer_to_objref(pointer(dp.ps_data))
            isa(ps_data,InsideSolidData)&&continue
            ap = Base.unsafe_wrap(
                Vector{Cdouble},
                mirror_slope_pointers[i],
                vs_num * DIM * NDF+DIM*(DIM+2),
            )
            get_mirror_slope_inner!(ps_data, ap)
        end
    end
end
function ghost_data_alloc(N::Int, ghost::P_pxest_ghost_t)
    Ptr{Cdouble}(sc_malloc(P4est.package_id(), sizeof(Cdouble) * N * PointerWrapper(ghost).ghosts.elem_count[]))
end
function ghost_data_alloc(::Type{T},N::Int, ghost::P_pxest_ghost_t) where{T}
    Ptr{T}(sc_malloc(P4est.package_id(), sizeof(Cdouble) * N * PointerWrapper(ghost).ghosts.elem_count[]))
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

"""
$(TYPEDSIGNATURES)
Update `df` in [`Ghost_VsData`](@ref).
"""
function data_exchange!(p4est::P_pxest_t, ka::KA{DIM,NDF}) where{DIM,NDF}
    MPI.Comm_size(MPI.COMM_WORLD) == 1 && return nothing
    update_mirror_data!(p4est, ka)
    amr_exchange!(
        p4est,
        ka.kinfo.forest.ghost,
        ka.kdata.ghost.ghost_pointers.ghost_datas,
        ka.kdata.ghost.ghost_pointers.mirror_data_pointers,
        size_Ghost_Data(ka.kinfo.status.max_vs_num,DIM,NDF),
    )
end

"""
$(TYPEDSIGNATURES)
Update `df` in [`Ghost_VsData`](@ref) for the update of immersed boundaries. 
    Currently, the communication occurs through all ghost layers. The difference between [`data_exchange!`](@ref) 
    is that only immersed-boundary-related ghost cells' `mirror_data` is updated.
"""
function solid_exchange!(p4est::P_pxest_t, ka::KA{DIM,NDF}) where{DIM,NDF}
    MPI.Comm_size(MPI.COMM_WORLD) == 1 && return nothing
    isempty(ka.kinfo.config.IB) && return nothing
    update_solid_mirror_data!(p4est, ka)
    amr_exchange!(
        p4est,
        ka.kinfo.forest.ghost,
        ka.kdata.ghost.ghost_pointers.ghost_datas,
        ka.kdata.ghost.ghost_pointers.mirror_data_pointers,
        size_Ghost_Data(ka.kinfo.status.max_vs_num,DIM,NDF),
    )
end

"""
$(TYPEDSIGNATURES)
Update `sdf` in [`Ghost_VsData`](@ref).
"""
function slope_exchange!(p4est::P_pxest_t, ka::KA{DIM,NDF}) where{DIM,NDF}
    MPI.Comm_size(MPI.COMM_WORLD) == 1 && return nothing
    update_mirror_slope!(p4est, ka)
    amr_exchange!(
        p4est,
        ka.kinfo.forest.ghost,
        ka.kdata.ghost.ghost_pointers.ghost_slopes,
        ka.kdata.ghost.ghost_pointers.mirror_slope_pointers,
        size_Ghost_Slope(ka.kinfo.status.max_vs_num,DIM,NDF),
    )
end
function initialize_ghost_pointers(p4est, kinfo::KInfo{DIM,NDF}) where{DIM,NDF}
    mirror_data_pointers = get_mirror_data(p4est, kinfo)
    ghost_datas = amr_exchange(
        p4est,
        kinfo.forest.ghost,
        mirror_data_pointers,
        size_Ghost_Data(kinfo.status.max_vs_num,DIM,NDF),
    )
    mirror_slope_pointers = get_mirror_slope(p4est, kinfo)
    ghost_slopes = amr_exchange(
        p4est,
        kinfo.forest.ghost,
        mirror_slope_pointers,
        size_Ghost_Slope(kinfo.status.max_vs_num,DIM,NDF),
    )
    mirror_structure_pointers = get_mirror_structure(p4est, kinfo)
    ghost_structures = amr_exchange(
        p4est,
        kinfo.forest.ghost,
        mirror_structure_pointers,
        size_Ghost_VS_Structure(kinfo.status.max_vs_num,DIM),
    )
    return GhostPointers(
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
function initialize_ghost_wrap(kinfo::KInfo{DIM,NDF}, ghost_pointers::GhostPointers) where{DIM,NDF}
    ghost_wrap =
        Array{AbstractGhostPsData}(undef, PointerWrapper(kinfo.forest.ghost).ghosts.elem_count[])
    global_vs_num = kinfo.status.max_vs_num
    for i in eachindex(ghost_wrap)
        pq = iPointerWrapper(PointerWrapper(kinfo.forest.ghost).ghosts, p4est_quadrant_t, i - 1)
        owner_rank = pq.p.piggy1.owner_rank[]
        which_tree = pq.p.piggy3.which_tree[];local_num = pq.p.piggy3.local_num[]
        quadid = which_tree*2^(DIM*kinfo.config.solver.AMR_PS_MAXLEVEL)+local_num
        p = ghost_data_ptr(
            size_Ghost_Data(global_vs_num,DIM,NDF),
            ghost_pointers.ghost_datas,
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
        bound_enc = Int(Base.unsafe_load(p + offset * sizeof(Cdouble)))
        if w[1]==Inf
            ghost_wrap[i] = GhostInsideSolidData{DIM,NDF}(bound_enc,midpoint,ds)
            continue
        end
        offset += 1
        # midpoint_vs = Base.unsafe_wrap(Matrix{Cdouble},Ptr{Cdouble}(pointer(ghost_data.micro)),(vs_num,DIM))
        df = Base.unsafe_wrap(Matrix{Cdouble}, p + offset * sizeof(Cdouble), (vs_num, NDF))
        p = ghost_data_ptr(
            size_Ghost_Slope(global_vs_num,DIM,NDF),
            ghost_pointers.ghost_slopes,
            i - 1,
        )
        sdf = Base.unsafe_wrap(Array{Cdouble}, p, (vs_num, NDF, DIM))
        sw = Base.unsafe_wrap(Matrix{Cdouble}, p + vs_num * NDF * DIM * sizeof(Cdouble), (DIM+2, DIM))
        p = ghost_data_ptr(
            size_Ghost_VS_Structure(global_vs_num,DIM),
            ghost_pointers.ghost_structures,
            i - 1,
        )
        offset = 0
        weight = Base.unsafe_wrap(Vector{Cdouble}, p, vs_num)
        offset += global_vs_num
        level = Base.unsafe_wrap(Vector{Cdouble}, p + offset * sizeof(Cdouble), vs_num)
        offset += global_vs_num
        midpoint_vs =
            Base.unsafe_wrap(Matrix{Cdouble}, p + offset * sizeof(Cdouble), (vs_num, DIM))
        vs_data = Ghost_VsData{DIM,NDF}(vs_num, Int8.(round.(level)), weight, midpoint_vs, df, sdf)
        ghost_wrap[i] = GhostPsData{DIM,NDF}(owner_rank,quadid,bound_enc,ds, midpoint, w, sw, vs_data)
    end
    return ghost_wrap
end

function update_ghost!(p4est::Ptr{p4est_t}, ka::KA)
    ghost_pointers = ka.kdata.ghost.ghost_pointers
    kinfo = ka.kinfo
    finalize_ghost!(ghost_pointers)
    p4est_ghost_destroy(kinfo.forest.ghost)
    kinfo.forest.ghost = p4est_ghost_new(p4est, P4EST_CONNECT_FULL)
    ka.kdata.ghost.ghost_pointers = initialize_ghost_pointers(p4est, kinfo)
    ka.kdata.ghost.ghost_wrap = initialize_ghost_wrap(kinfo, ka.kdata.ghost.ghost_pointers)
end
function update_ghost!(p4est::Ptr{p8est_t}, ka::KA)
    ghost_pointers = ka.kdata.ghost.ghost_pointers
    kinfo = ka.kinfo
    finalize_ghost!(ghost_pointers)
    p8est_ghost_destroy(kinfo.forest.ghost)
    kinfo.forest.ghost = p8est_ghost_new(p4est, P8EST_CONNECT_FULL)
    ka.kdata.ghost.ghost_pointers = initialize_ghost_pointers(p4est, kinfo)
    ka.kdata.ghost.ghost_wrap = initialize_ghost_wrap(kinfo, ka.kdata.ghost.ghost_pointers)
end
