# Inner
"""
$(TYPEDSIGNATURES)
"""
function update_criterion!(
    ::Val{1},
    ::Val{1},
    ps_data::PsData,
    kinfo::KInfo{DIM,NDF},
    Ldata::AbstractVector,
    Rdata::AbstractVector,
    dir::Integer,
    ws_swL::Vector{Float64},
    ws_swR::Vector{Float64},
) where {DIM,NDF}
    ds = ps_data.ds[dir]
    update_Lohner_inner_ps!(ps_data, Ldata, Rdata, ds, ds, dir, ws_swL, ws_swR, kinfo)
    return nothing
end
"""
$(TYPEDSIGNATURES)
"""
function update_criterion!(
    ::Val{1},
    ::NeighborNum,
    ps_data::PsData,
    kinfo::KInfo{DIM,NDF},
    Ldata::AbstractVector,
    Rdata::AbstractVector,
    dir::Integer,
    ws_swL::Vector{Float64},
    ws_swR::Vector{Float64},
) where {DIM,NDF}
    ds = ps_data.ds[dir]
    update_Lohner_inner_ps!(
        ps_data,
        Ldata,
        Rdata,
        ds,
        0.75 * ds,
        dir,
        ws_swL,
        ws_swR,
        kinfo,
    )
end
"""
$(TYPEDSIGNATURES)
"""
function update_criterion!(
    ::NeighborNum,
    ::Val{1},
    ps_data::PsData,
    kinfo::KInfo{DIM,NDF},
    Ldata::AbstractVector,
    Rdata::AbstractVector,
    dir::Integer,
    ws_swL::Vector{Float64},
    ws_swR::Vector{Float64},
) where {DIM,NDF}
    ds = ps_data.ds[dir]
    update_Lohner_inner_ps!(
        ps_data,
        Ldata,
        Rdata,
        0.75 * ds,
        ds,
        dir,
        ws_swL,
        ws_swR,
        kinfo,
    )
end
"""
$(TYPEDSIGNATURES)
"""
function update_criterion!(
    ::NeighborNum,
    ::NeighborNum,
    ps_data::PsData,
    kinfo::KInfo{DIM,NDF},
    Ldata::AbstractVector,
    Rdata::AbstractVector,
    dir::Integer,
    ws_swL::Vector{Float64},
    ws_swR::Vector{Float64},
) where {DIM,NDF}
    ds = ps_data.ds[dir]
    update_Lohner_inner_ps!(
        ps_data,
        Ldata,
        Rdata,
        0.75 * ds,
        0.75 * ds,
        dir,
        ws_swL,
        ws_swR,
        kinfo,
    )
end
"""
$(TYPEDSIGNATURES)
"""
function update_criterion!(
    ::Val{-1},
    ::Val{-1},
    ps_data::PsData,
    kinfo::KInfo{DIM,NDF},
    Ldata::AbstractVector,
    Rdata::AbstractVector,
    dir::Integer,
    ws_sL::Matrix{Float64},
    ws_sR::Matrix{Float64},
    ws_swL::Vector{Float64},
    ws_swR::Vector{Float64},
) where {DIM,NDF}
    ds = ps_data.ds[dir]
    update_Lohner_inner_ps!(
        ps_data,
        Ldata,
        Rdata,
        1.5 * ds,
        1.5 * ds,
        dir,
        ws_swL,
        ws_swR,
        kinfo,
    )
end
"""
$(TYPEDSIGNATURES)
"""
function update_criterion!(
    ::Val{1},
    ::Val{-1},
    ps_data::PsData,
    kinfo::KInfo{DIM,NDF},
    Ldata::AbstractVector,
    Rdata::AbstractVector,
    dir::Integer,
    ws_swL::Vector{Float64},
    ws_swR::Vector{Float64},
) where {DIM,NDF}
    ds = ps_data.ds[dir]
    update_Lohner_inner_ps!(ps_data, Ldata, Rdata, ds, 1.5 * ds, dir, ws_swL, ws_swR, kinfo)
end
"""
$(TYPEDSIGNATURES)
"""
function update_criterion!(
    ::Val{-1},
    ::Val{1},
    ps_data::PsData,
    kinfo::KInfo{DIM,NDF},
    Ldata::AbstractVector,
    Rdata::AbstractVector,
    dir::Integer,
    ws_swL::Vector{Float64},
    ws_swR::Vector{Float64},
) where {DIM,NDF}
    ds = ps_data.ds[dir]
    update_Lohner_inner_ps!(ps_data, Ldata, Rdata, 1.5 * ds, ds, dir, ws_swL, ws_swR, kinfo)
end
"""
$(TYPEDSIGNATURES)
"""
function update_criterion!(
    ::NeighborNum,
    ::Val{-1},
    ps_data::PsData,
    kinfo::KInfo{DIM,NDF},
    Ldata::AbstractVector,
    Rdata::AbstractVector,
    dir::Integer,
    ws_swL::Vector{Float64},
    ws_swR::Vector{Float64},
) where {DIM,NDF}
    ds = ps_data.ds[dir]
    update_Lohner_inner_ps!(
        ps_data,
        Ldata,
        Rdata,
        0.75 * ds,
        1.5 * ds,
        dir,
        ws_swL,
        ws_swR,
        kinfo,
    )
end
"""
$(TYPEDSIGNATURES)
"""
function update_criterion!(
    ::Val{-1},
    ::NeighborNum,
    ps_data::PsData,
    kinfo::KInfo{DIM,NDF},
    Ldata::AbstractVector,
    Rdata::AbstractVector,
    dir::Integer,
    ws_swL::Vector{Float64},
    ws_swR::Vector{Float64},
) where {DIM,NDF}
    ds = ps_data.ds[dir]
    update_Lohner_inner_ps!(
        ps_data,
        Ldata,
        Rdata,
        1.5 * ds,
        0.75 * ds,
        dir,
        ws_swL,
        ws_swR,
        kinfo,
    )
end

# Boundary
"""
$(TYPEDSIGNATURES)
"""
function update_criterion!(
    ::Val{0},
    ::Val{1},
    ps_data::PsData,
    kinfo::KInfo{DIM,NDF},
    Ldata::AbstractVector,
    Rdata::AbstractVector,
    dir::Integer,
    ws_swL::Vector{Float64},
    ws_swR::Vector{Float64},
) where {DIM,NDF}
    ps_data.lohner[:, dir] .= 0.0
end
"""
$(TYPEDSIGNATURES)
"""
function update_criterion!(
    ::Val{1},
    ::Val{0},
    ps_data::PsData,
    kinfo::KInfo{DIM,NDF},
    Ldata::AbstractVector,
    Rdata::AbstractVector,
    dir::Integer,
    ws_swL::Vector{Float64},
    ws_swR::Vector{Float64},
) where {DIM,NDF}
    ps_data.lohner[:, dir] .= 0.0
end
"""
$(TYPEDSIGNATURES)
"""
function update_criterion!(
    ::NeighborNum,
    ::Val{0},
    ps_data::PsData,
    kinfo::KInfo{DIM,NDF},
    Ldata::AbstractVector,
    Rdata::AbstractVector,
    dir::Integer,
    ws_swL::Vector{Float64},
    ws_swR::Vector{Float64},
) where {DIM,NDF}
    ps_data.lohner[:, dir] .= 0.0
end
"""
$(TYPEDSIGNATURES)
"""
function update_criterion!(
    ::Val{0},
    ::NeighborNum,
    ps_data::PsData,
    kinfo::KInfo{DIM,NDF},
    Ldata::AbstractVector,
    Rdata::AbstractVector,
    dir::Integer,
    ws_swL::Vector{Float64},
    ws_swR::Vector{Float64},
) where {DIM,NDF}
    ps_data.lohner[:, dir] .= 0.0
end
"""
$(TYPEDSIGNATURES)
"""
function update_criterion!(
    ::Val{0},
    ::Val{-1},
    ps_data::PsData,
    kinfo::KInfo{DIM,NDF},
    Ldata::AbstractVector,
    Rdata::AbstractVector,
    dir::Integer,
    ws_swL::Vector{Float64},
    ws_swR::Vector{Float64},
) where {DIM,NDF}
    ps_data.lohner[:, dir] .= 0.0
end
"""
$(TYPEDSIGNATURES)
"""
function update_criterion!(
    ::Val{-1},
    ::Val{0},
    ps_data::PsData,
    kinfo::KInfo{DIM,NDF},
    Ldata::AbstractVector,
    Rdata::AbstractVector,
    dir::Integer,
    ws_swL::Vector{Float64},
    ws_swR::Vector{Float64},
) where {DIM,NDF}
    ps_data.lohner[:, dir] .= 0.0
end
function update_criterion!(ka::KA{DIM,NDF}) where {DIM,NDF}
    trees = ka.kdata.field.trees
    kinfo = ka.kinfo
    ws_swL = zeros(Float64, DIM + 2)
    ws_swR = zeros(Float64, DIM + 2)

    @inbounds for i in eachindex(trees.data)
        @inbounds for j in eachindex(trees.data[i])
            ps_data = trees.data[i][j]
            isa(ps_data, InsideSolidData) && continue
            ps_data.bound_enc<0 && continue # solid_cells
            neighbor = ps_data.neighbor
            for dir = 1:DIM
                iL = 2 * dir - 1
                iR = 2 * dir
                update_criterion!(
                    Val(neighbor.state[iL]),
                    Val(neighbor.state[iR]),
                    ps_data,
                    kinfo,
                    neighbor.data[iL],
                    neighbor.data[iR],
                    dir,
                    ws_swL,
                    ws_swR,
                )
            end
        end
    end
    apply_amr_buffer!(ka)
    return nothing
end

"""
$(TYPEDSIGNATURES)
One-cell buffer layer: any cell whose neighbor's Löhner sensor exceeds the
refine threshold has its own `lohner` inflated above the threshold, so it
will also be flagged for refinement. Decide-then-apply keeps the buffer
strictly one layer thick. Ghost neighbours' flags are synced via
[`lohner_flag_exchange!`](@ref).
"""
function apply_amr_buffer!(ka::KA{DIM,NDF}) where {DIM,NDF}
    trees = ka.kdata.field.trees
    threshold = ka.kinfo.config.solver.ADAPT_COEFFI_PS
    inflate = 2.0 * threshold

    ghost_flags = lohner_flag_exchange!(ka.kinfo.forest.p4est, ka)

    buffer = [falses(length(t)) for t in trees.data]
    @inbounds for i in eachindex(trees.data)
        for j in eachindex(trees.data[i])
            ps_data = trees.data[i][j]
            isa(ps_data, InsideSolidData) && continue
            ps_data.bound_enc < 0 && continue
            ps_sensor(ps_data) > threshold && continue
            neighbor = ps_data.neighbor
            flagged = false
            for k in eachindex(neighbor.data)
                neighbor.state[k] == 0 && continue
                for nb in neighbor.data[k]
                    if isa(nb, PsData)
                        nb.bound_enc < 0 && continue
                        if ps_sensor(nb) > threshold
                            flagged = true
                            break
                        end
                    elseif isa(nb, GhostPsData)
                        nb.bound_enc < 0 && continue
                        if get(ghost_flags, objectid(nb), false)
                            flagged = true
                            break
                        end
                    end
                end
                flagged && break
            end
            buffer[i][j] = flagged
        end
    end

    @inbounds for i in eachindex(trees.data)
        for j in eachindex(trees.data[i])
            buffer[i][j] && (trees.data[i][j].lohner .= inflate)
        end
    end
    return nothing
end

function vs_merge!(
    sdf::AbstractArray,
    sdf_n::AbstractArray,
    level::Vector,
    level_n::Vector,
    ::KA{DIM},
) where {DIM}
    j = 1
    flag = 0.0
    for i in axes(sdf, 1)
        if level[i] == level_n[j]
            @. sdf[i, :, :] += @views sdf_n[j, :, :]
            j += 1
        elseif level[i] < level_n[j]
            while flag != 1.0
                @. sdf[i, :, :] += @views sdf_n[j, :, :] / 2^(DIM * (level_n[j] - level[i]))
                flag += 1 / 2^(DIM * (level_n[j] - level[i]))
                j += 1
            end
            flag = 0.0
        else
            @. sdf[i, :, :] += @views sdf_n[j, :, :]
            flag += 1 / 2^(DIM * (level[i] - level_n[j]))
            if flag == 1.0
                j += 1
                flag = 0.0
            end
        end
    end
end
function vs_merge!(
    sdf::AbstractMatrix,
    sdf_n::AbstractMatrix,
    level::Vector,
    level_n::Vector,
    ::KA{DIM},
) where {DIM}
    j = 1
    flag = 0.0
    for i in axes(sdf, 1)
        if level[i] == level_n[j]
            @. sdf[i, :] += @views sdf_n[j, :]
            j += 1
        elseif level[i] < level_n[j]
            while flag != 1.0
                @. sdf[i, :] += @views sdf_n[j, :] / 2^(DIM * (level_n[j] - level[i]))
                flag += 1 / 2^(DIM * (level_n[j] - level[i]))
                j += 1
            end
            flag = 0.0
        else
            @. sdf[i, :] += @views sdf_n[j, :]
            flag += 1 / 2^(DIM * (level[i] - level_n[j]))
            if flag == 1.0
                j += 1
                flag = 0.0
            end
        end
    end
end


function user_defined_ps_refine_flag(forest::P_pxest_t, which_tree, quadrant)
    GC.@preserve forest which_tree quadrant begin
        fp = PointerWrapper(forest)
        qp = PointerWrapper(quadrant)
        kinfo = unsafe_pointer_to_objref(pointer(fp.user_pointer))
        refine_fn = kinfo.config.user_defined.static_ps_refine_flag
        ds, midpoint = quad_to_cell(fp, which_tree, qp)
        if refine_fn(midpoint, ds, kinfo, qp.level[])
            return Cint(1)
        else
            return Cint(0)
        end
    end
end
function user_defined_ps_refine!(p4est::Ptr{p4est_t}, kinfo::KInfo)
    p4est_refine_ext(
        p4est,
        1,
        kinfo.config.solver.AMR_PS_MAXLEVEL,
        @cfunction(
            user_defined_ps_refine_flag,
            Cint,
            (Ptr{p4est_t}, p4est_topidx_t, Ptr{p4est_quadrant_t})
        ),
        C_NULL,
        C_NULL,
    )
end
function user_defined_ps_refine!(p4est::Ptr{p8est_t}, kinfo::KInfo)
    p8est_refine_ext(
        p4est,
        1,
        kinfo.config.solver.AMR_PS_MAXLEVEL,
        @cfunction(
            user_defined_ps_refine_flag,
            Cint,
            (Ptr{p8est_t}, p4est_topidx_t, Ptr{p8est_quadrant_t})
        ),
        C_NULL,
        C_NULL,
    )
end

function initialize_MeshData!(p4est::P_pxest_t, kinfo::KInfo)
    fp = PointerWrapper(p4est)
    trees = [MeshData[] for _ = 1:(fp.last_local_tree[]-fp.first_local_tree[]+1)]
    p_data = pointer_from_objref(trees)
    GC.@preserve trees AMR_volume_iterate(p4est; user_data = p_data) do ip, data, dp
        treeid = ip.treeid[]-ip.p4est.first_local_tree[]+1
        mesh_data = MeshData()
        dp[] = P4estPsData(pointer_from_objref(mesh_data))
        trees = unsafe_pointer_to_objref(data)
        push!(trees[treeid], mesh_data)
    end
    return trees
end

function search_radius_flag!(forest::P_pxest_t, which_tree, quadrant)
    GC.@preserve forest which_tree quadrant begin
        fp = PointerWrapper(forest)
        kinfo, _ = unsafe_pointer_to_objref(pointer(fp.user_pointer))
        qp = PointerWrapper(quadrant)
        dp = PointerWrapper(P4estPsData, qp.p.user_data[])
        mesh_data = unsafe_pointer_to_objref(pointer(dp.ps_data))
        ibs = kinfo.config.IB
        ds, midpoint = quad_to_cell(fp, which_tree, qp)
        for i in eachindex(ibs)
            search_radius_flag!(i, ibs[i], midpoint, ds, mesh_data) && return Cint(1)
        end
        return Cint(0)
    end
end
function search_radius_refine_replace(
    forest::P_pxest_t,
    which_tree,
    num_out,
    out_quads::Ptr{T},
    num_in,
    in_quads,
) where {T}
    GC.@preserve forest which_tree num_out out_quads num_in in_quads begin
        fp = PointerWrapper(forest)
        kinfo, trees = unsafe_pointer_to_objref(pointer(fp.user_pointer))
        out_quads_wrap = unsafe_wrap(Vector{T}, out_quads, num_out)
        in_quads_wrap = unsafe_wrap(Vector{T}, in_quads, num_in)
        search_radius_refine_replace(
            forest,
            out_quads_wrap,
            in_quads_wrap,
            which_tree,
            trees,
            kinfo,
        )
        return nothing
    end
end
function search_radius_refine!(p4est::Ptr{p4est_t}, kinfo::KInfo)
    p4est_refine_ext(
        p4est,
        1,
        kinfo.config.solver.AMR_PS_MAXLEVEL,
        @cfunction(
            search_radius_flag!,
            Cint,
            (Ptr{p4est_t}, p4est_topidx_t, Ptr{p4est_quadrant_t})
        ),
        C_NULL,
        @cfunction(
            search_radius_refine_replace,
            Cvoid,
            (
                Ptr{p4est_t},
                p4est_topidx_t,
                Cint,
                Ptr{Ptr{p4est_quadrant_t}},
                Cint,
                Ptr{Ptr{p4est_quadrant_t}},
            )
        )
    )
end
function search_radius_refine!(p4est::Ptr{p8est_t}, kinfo::KInfo)
    p8est_refine_ext(
        p4est,
        1,
        kinfo.config.solver.AMR_PS_MAXLEVEL,
        @cfunction(
            search_radius_flag!,
            Cint,
            (Ptr{p8est_t}, p4est_topidx_t, Ptr{p8est_quadrant_t})
        ),
        C_NULL,
        @cfunction(
            search_radius_refine_replace,
            Cvoid,
            (
                Ptr{p8est_t},
                p4est_topidx_t,
                Cint,
                Ptr{Ptr{p8est_quadrant_t}},
                Cint,
                Ptr{Ptr{p8est_quadrant_t}},
            )
        )
    )
end

function search_radius_refine_replace(p4est, out_quads, in_quads, which_tree, trees, kinfo)
    fp = PointerWrapper(p4est)
    out_qp = PointerWrapper(out_quads[1])
    dp = PointerWrapper(P4estPsData, out_qp.p.user_data[])
    mesh_data = unsafe_pointer_to_objref(pointer(dp.ps_data))
    ib = kinfo.config.IB[mesh_data.in_search_radius]
    treeid = which_tree - fp.first_local_tree[]+1
    datas = trees[treeid]
    index = findfirst(x->x===mesh_data, datas)
    deleteat!(datas, index)
    if out_qp.level[]==kinfo.config.solver.AMR_PS_MAXLEVEL-1
        for i in eachindex(in_quads)
            qp = PointerWrapper(in_quads[i])
            dp = PointerWrapper(P4estPsData, qp.p.user_data[])
            mesh_data_i = deepcopy(mesh_data)
            dp[] = P4estPsData(pointer_from_objref(mesh_data_i))
            insert!(datas, index - 1 + i, mesh_data_i)
        end
    else
        for i in eachindex(in_quads)
            qp = PointerWrapper(in_quads[i])
            dp = PointerWrapper(P4estPsData, qp.p.user_data[])
            mesh_data_i = MeshData()
            dp[] = P4estPsData(pointer_from_objref(mesh_data_i))
            insert!(datas, index - 1 + i, mesh_data_i)
        end
    end
end
# function boundary_flag(::Domain,::AbstractVector,::AbstractVector,::KInfo) # Domain boundary flag
#     return false
# end
function boundary_flag(
    boundary::Domain,
    midpoint::AbstractVector,
    ds::AbstractVector,
    kinfo::KInfo,
) # Domain boundary flag
    !boundary.refine && return false
    index = Int(floor((boundary.id-1)/2))+1
    abs(midpoint[index] - kinfo.config.geometry[boundary.id]) < ds[index] && return true
end
function domain_flag(kinfo::KInfo, midpoint::AbstractVector, ds::AbstractVector) # domain flag for physical space dynamic adaptive
    domains = kinfo.config.domain
    for domain in domains
        boundary_flag(domain, midpoint, ds, kinfo) && return true
    end
    return false
end
function boundary_flag(midpoint::AbstractVector, ds::AbstractVector, kinfo::KInfo) # Is the cell at midpoint 
    domain = kinfo.config.domain
    boundary = kinfo.config.IB
    for i in eachindex(domain)
        boundary_flag(domain[i], midpoint, ds, kinfo) && return true
    end
    return solid_cell_flag(boundary, midpoint, ds, kinfo)
end
function ps_copy(data::PsData{DIM,NDF}) where {DIM,NDF}
    p = PsData(
        DIM,
        NDF;
        ds = copy(data.ds),
        midpoint = copy(data.midpoint),
        w = copy(data.w),
        sw = copy(data.sw),
        lohner = copy(data.lohner),
        vs_data = deepcopy(data.vs_data),
        prim = copy(data.prim),
        flux = copy(data.flux),
        # neighbor = data.neighbor
    )
    return p
end

function ps_merge(Odatas::Vector, index::Int, kinfo::KInfo{DIM,NDF}) where {DIM,NDF}
    data = Odatas[index]
    vs_data = data.vs_data;
    vs_num = vs_data.vs_num
    w_new = [sum([x.w[i] for x in Odatas]) for i = 1:(DIM+2)] ./ 2^DIM
    p = PsData(
        DIM,
        NDF;
        ds = data.ds * 2.0,
        midpoint = data.midpoint - 0.5 * data.ds .* RMT[DIM][index],
        w = w_new,
        prim = get_prim(w_new, kinfo),
        sw = [sum([x.sw[i, j] for x in Odatas]) for i = 1:(DIM+2), j = 1:DIM] ./ 2^DIM,
        flux = [sum([x.flux[i] for x in Odatas]) for i = 1:(DIM+2)] ./ 2^DIM,
        vs_data = VsData{DIM,NDF}(
            vs_num,
            copy(vs_data.level),
            copy(vs_data.weight),
            copy(vs_data.midpoint),
            zeros(vs_num, NDF),
            zeros(vs_num, NDF, DIM),
            zeros(vs_num, NDF),
        ),
    )
    p.neighbor.state[1] = BALANCE_FLAG
    p.neighbor.data[1] = Odatas
    return p
end
function ps_refine_flag(::InsideSolidData, ::Int8, ::KA)
    Cint(0)
end
function ps_refine_flag(forest::P_pxest_t, which_tree, quadrant)
    GC.@preserve forest which_tree quadrant begin
        fp = PointerWrapper(forest)
        qp = PointerWrapper(quadrant)
        level = qp.level[]
        dp = PointerWrapper(P4estPsData, qp.p.user_data[])
        ps_data = unsafe_pointer_to_objref(pointer(dp.ps_data))
        ka = unsafe_pointer_to_objref(pointer(fp.user_pointer))
        ps_refine_flag(ps_data, level, ka)
    end
end

"""
$(TYPEDSIGNATURES)
Replacement function in a refinement process.
"""
function ps_replace!(::Val{1}, out_quad, in_quads, which_tree, ka::KA{DIM}) where {DIM}# refine replace
    trees = ka.kdata.field.trees
    treeid = Int(which_tree) - trees.offset
    datas = trees.data[treeid]
    pw_out_quad = PointerWrapper(out_quad[1])
    Odata = unsafe_pointer_to_objref(
        pointer(PointerWrapper(P4estPsData, pw_out_quad.p.user_data[]).ps_data),
    )
    index = findfirst(x -> x === Odata, datas)
    deleteat!(datas, index)
    flag = Odata.neighbor.state[1]==BALANCE_FLAG
    for i = 1:(2^DIM)
        pw_in_quad = PointerWrapper(in_quads[i])
        dp = PointerWrapper(P4estPsData, pw_in_quad.p.user_data[])
        if flag
            ps_data = Odata.neighbor.data[1][i]
        else
            ps_data = ps_copy(Odata)
            ps_data.ds .*= 0.5
            vs_data = ps_data.vs_data
            @. ps_data.midpoint += 0.5 * ps_data.ds * RMT[DIM][i]
            for k in axes(vs_data.df, 2)
                for j in axes(vs_data.df, 1)
                    @views vs_data.df[j, k]+=min(
                        abs(
                            vs_data.df[j, k]/(
                                dot(ps_data.ds, abs.(vs_data.sdf[j, k, :]))+EPS
                            ),
                        ),
                        1.0,
                    )*dot(vs_data.sdf[j, k, :], 0.5 * ps_data.ds .* RMT[DIM][i])
                end
            end
            ps_data.w = calc_w0(ps_data)
            ps_data.prim = get_prim(ps_data.w, ka.kinfo)
        end
        insert!(datas, index - 1 + i, ps_data)
        dp[] = P4estPsData(pointer_from_objref(ps_data))
    end
    return nothing
end
"""
$(TYPEDSIGNATURES)
Replacement function in a coarsening process.
"""
function ps_replace!(
    ::ChildNum,
    out_quad,
    in_quads,
    which_tree,
    ka::KA{DIM,NDF},
) where {DIM,NDF} # coarsen replace, average or interpolate? Currently interpolation strategy is adopted. If my memory serves me right, problems came out with average most likely due to the iterative balance process.
    trees = ka.kdata.field.trees
    treeid = Int(which_tree) - trees.offset
    datas = trees.data[treeid]
    pw_in_quad = PointerWrapper(in_quads[1])
    dp = PointerWrapper(P4estPsData, pw_in_quad.p.user_data[])
    Odatas = Vector{AbstractPsData{DIM,NDF}}(undef, 2^DIM)
    for i = 1:(2^DIM)
        pw_out_quad = PointerWrapper(out_quad[i])
        Odatas[i] = unsafe_pointer_to_objref(
            pointer(PointerWrapper(P4estPsData, pw_out_quad.p.user_data[]).ps_data),
        )
    end
    if any(x->x<0, [x.bound_enc for x in Odatas])
        @error `solid_cells coarsen!`
    end
    _, index = findmin([x.vs_data.vs_num for x in Odatas])
    ps_data = ps_merge(Odatas, index, ka.kinfo)
    vs_data = ps_data.vs_data
    for i in eachindex(Odatas)
        vs_data_n = Odatas[i].vs_data
        vs_merge!(vs_data.sdf, vs_data_n.sdf, vs_data.level, vs_data_n.level, ka)
        vs_merge!(vs_data.df, vs_data_n.df, vs_data.level, vs_data_n.level, ka)
    end
    vs_data.df./=2^DIM;
    vs_data.sdf./=2^DIM
    ps_data.neighbor.state[1] = BALANCE_FLAG
    ps_data.neighbor.data[1] = Odatas
    index = findfirst(x -> x === Odatas[1], datas)
    deleteat!(datas, index:(index+2^DIM-1))
    insert!(datas, index, ps_data)
    dp[] = P4estPsData(pointer_from_objref(ps_data))
    return nothing
end
function p4est_replace(
    forest::T1,
    which_tree,
    num_out,
    out_quads::Ptr{T2},
    num_in,
    in_quads,
) where {T1<:P_pxest_t,T2<:P_pxest_quadrant_t}
    GC.@preserve forest which_tree num_out out_quads num_in in_quads begin
        fp = PointerWrapper(forest)
        ka = unsafe_pointer_to_objref(pointer(fp.user_pointer))
        out_quads_wrap = unsafe_wrap(Vector{T2}, out_quads, num_out)
        in_quads_wrap = unsafe_wrap(Vector{T2}, in_quads, num_in)
        ps_replace!(Val(Int(num_out)), out_quads_wrap, in_quads_wrap, which_tree, ka)
        return nothing
    end
end
function pre_ps_replace(
    forest::T1,
    which_tree,
    num_out,
    out_quads::Ptr{T2},
    num_in,
    in_quads,
) where {T1<:P_pxest_t,T2<:P_pxest_quadrant_t}
    GC.@preserve forest which_tree num_out out_quads num_in in_quads begin
        fp = PointerWrapper(forest)
        _, trees = unsafe_pointer_to_objref(pointer(fp.user_pointer))
        out_quads_wrap = unsafe_wrap(Vector{T2}, out_quads, num_out)
        in_quads_wrap = unsafe_wrap(Vector{T2}, in_quads, num_in)
        pre_ps_replace(
            Val(Int(num_out)),
            forest,
            out_quads_wrap,
            in_quads_wrap,
            which_tree,
            trees,
        )
        return nothing
    end
end
function pre_ps_replace(::Val{1}, p4est, out_quad, in_quads, which_tree, trees)# refine replace
    fp = PointerWrapper(p4est)
    treeid = Int(which_tree) - fp.first_local_tree[]+1
    datas = trees[treeid]
    pw_out_quad = PointerWrapper(out_quad[1])
    Odata = unsafe_pointer_to_objref(
        pointer(PointerWrapper(P4estPsData, pw_out_quad.p.user_data[]).ps_data),
    )
    index = findfirst(x -> x === Odata, datas)
    deleteat!(datas, index)
    for i in eachindex(in_quads)
        pw_in_quad = PointerWrapper(in_quads[i])
        dp = PointerWrapper(P4estPsData, pw_in_quad.p.user_data[])
        mesh_data = deepcopy(Odata)
        insert!(datas, index - 1 + i, mesh_data)
        dp[] = P4estPsData(pointer_from_objref(mesh_data))
    end
end
function pre_ps_replace(
    ::ChildNum,
    p4est,
    out_quad,
    in_quads,
    which_tree,
    trees::Vector{Vector{T}},
) where {T} # coarsen replace, average or interpolate? Currently interpolation strategy is adopted. If my memory serves me right, problems came out with average most likely due to the iterative balance process.
    fp = PointerWrapper(p4est)
    treeid = Int(which_tree) - fp.first_local_tree[]+1
    datas = trees[treeid]
    pw_in_quad = PointerWrapper(in_quads[1])
    dp = PointerWrapper(P4estPsData, pw_in_quad.p.user_data[])
    Odatas = Vector{T}(undef, length(out_quad))
    for i in eachindex(out_quad)
        pw_out_quad = PointerWrapper(out_quad[i])
        Odatas[i] = unsafe_pointer_to_objref(
            pointer(PointerWrapper(P4estPsData, pw_out_quad.p.user_data[]).ps_data),
        )
    end
    mesh_data = deepcopy(Odatas[1])
    index = findfirst(x -> x === Odatas[1], datas)
    deleteat!(datas, index:(index+length(out_quad)-1))
    insert!(datas, index, mesh_data)
    dp[] = P4estPsData(pointer_from_objref(mesh_data))
end

"""
$(TYPEDSIGNATURES)
"""
function ps_refine!(p4est::Ptr{p4est_t}, ka::KA; recursive = 0)
    p4est_refine_ext(
        p4est,
        recursive,
        ka.kinfo.config.solver.AMR_DYNAMIC_PS_MAXLEVEL,
        @cfunction(
            ps_refine_flag,
            Cint,
            (Ptr{p4est_t}, p4est_topidx_t, Ptr{p4est_quadrant_t})
        ),
        C_NULL,
        @cfunction(
            p4est_replace,
            Cvoid,
            (
                Ptr{p4est_t},
                p4est_topidx_t,
                Cint,
                Ptr{Ptr{p4est_quadrant_t}},
                Cint,
                Ptr{Ptr{p4est_quadrant_t}},
            )
        )
    )
end
"""
$(TYPEDSIGNATURES)
"""
function ps_refine!(p4est::Ptr{p8est_t}, ka::KA; recursive = 0)
    p8est_refine_ext(
        p4est,
        recursive,
        ka.kinfo.config.solver.AMR_PS_MAXLEVEL,
        @cfunction(
            ps_refine_flag,
            Cint,
            (Ptr{p8est_t}, p4est_topidx_t, Ptr{p8est_quadrant_t})
        ),
        C_NULL,
        @cfunction(
            p4est_replace,
            Cvoid,
            (
                Ptr{p8est_t},
                p4est_topidx_t,
                Cint,
                Ptr{Ptr{p8est_quadrant_t}},
                Cint,
                Ptr{Ptr{p8est_quadrant_t}},
            )
        )
    )
end


function pre_ps_coarsen_flag(forest::Ptr{p4est_t}, which_tree, quadrants)
    GC.@preserve forest which_tree quadrants begin
        DIM = 2
        fp = PointerWrapper(forest)
        quadrants_wrap = unsafe_wrap(Vector{Ptr{p4est_quadrant_t}}, quadrants, 2^DIM)
        kinfo, _ = unsafe_pointer_to_objref(pointer(fp.user_pointer))
        for i = 1:(2^DIM)
            quadrants_wrap[i]==C_NULL&&break
            qp = PointerWrapper(quadrants_wrap[i])
            dp = PointerWrapper(P4estPsData, qp.p.user_data[])
            mesh_data = unsafe_pointer_to_objref(pointer(dp.ps_data))
            ds, midpoint = quad_to_cell(fp, which_tree, qp)
            (
                kinfo.config.user_defined.static_ps_refine_flag(
                    midpoint,
                    ds,
                    kinfo,
                    qp.level[]-1,
                )&&!mesh_data.in_solid ||
                mesh_data.is_ghost_cell||(
                    mesh_data.in_search_radius!=0&&!mesh_data.in_solid
                )
            )&&return Cint(0)
        end
        return Cint(1)
    end
end
function pre_ps_coarsen_flag(forest::Ptr{p8est_t}, which_tree, quadrants)
    GC.@preserve forest which_tree quadrants begin
        DIM = 3
        fp = PointerWrapper(forest)
        quadrants_wrap = unsafe_wrap(Vector{Ptr{p8est_quadrant_t}}, quadrants, 2^DIM)
        kinfo, _ = unsafe_pointer_to_objref(pointer(fp.user_pointer))
        for i = 1:(2^DIM)
            qp = PointerWrapper(quadrants_wrap[i])
            dp = PointerWrapper(P4estPsData, qp.p.user_data[])
            mesh_data = unsafe_pointer_to_objref(pointer(dp.ps_data))
            ds, midpoint = quad_to_cell(fp, which_tree, qp)
            (
                kinfo.config.user_defined.static_ps_refine_flag(
                    midpoint,
                    ds,
                    kinfo,
                    qp.level[]-1,
                )&&!mesh_data.in_solid ||
                mesh_data.is_ghost_cell||(
                    mesh_data.in_search_radius!=0&&!mesh_data.in_solid
                )
            )&&return Cint(0)
        end
        return Cint(1)
    end
end
function pre_ps_coarsen!(p4est::Ptr{p4est_t}; recursive = 0)
    p4est_coarsen_ext(
        p4est,
        recursive,
        1,
        @cfunction(
            pre_ps_coarsen_flag,
            Cint,
            (Ptr{p4est_t}, p4est_topidx_t, Ptr{Ptr{p4est_quadrant_t}})
        ),
        C_NULL,
        @cfunction(
            pre_ps_replace,
            Cvoid,
            (
                Ptr{p4est_t},
                p4est_topidx_t,
                Cint,
                Ptr{Ptr{p4est_quadrant_t}},
                Cint,
                Ptr{Ptr{p4est_quadrant_t}},
            )
        )
    )
end
function pre_ps_coarsen!(p4est::Ptr{p8est_t}; recursive = 0)
    p8est_coarsen_ext(
        p4est,
        recursive,
        0,
        @cfunction(
            pre_ps_coarsen_flag,
            Cint,
            (Ptr{p8est_t}, p4est_topidx_t, Ptr{Ptr{p8est_quadrant_t}})
        ),
        C_NULL,
        @cfunction(
            pre_ps_replace,
            Cvoid,
            (
                Ptr{p8est_t},
                p4est_topidx_t,
                Cint,
                Ptr{Ptr{p8est_quadrant_t}},
                Cint,
                Ptr{Ptr{p8est_quadrant_t}},
            )
        )
    )
end

function ps_coarsen_flag(forest::Ptr{p4est_t}, which_tree, quadrants)
    GC.@preserve forest which_tree quadrants begin
        DIM = 2
        fp = PointerWrapper(forest)
        quadrants_wrap = unsafe_wrap(Vector{Ptr{p4est_quadrant_t}}, quadrants, 2^DIM)
        levels = Vector{Int}(undef, 2^DIM)
        ps_datas = Vector{PsData}(undef, 2^DIM)
        for i = 1:(2^DIM)
            qp = PointerWrapper(quadrants_wrap[i])
            dp = PointerWrapper(P4estPsData, qp.p.user_data[])
            ps_data = unsafe_pointer_to_objref(pointer(dp.ps_data))
            isa(ps_data, InsideSolidData)&&return Cint(0)
            levels[i] = qp.level[]
            ps_datas[i] = ps_data
        end
        ka = unsafe_pointer_to_objref(pointer(fp.user_pointer))
        ps_coarsen_flag(ps_datas, levels, ka)
    end
end
function ps_coarsen_flag(forest::Ptr{p8est_t}, which_tree, quadrants)
    GC.@preserve forest which_tree quadrants begin
        DIM = 3
        fp = PointerWrapper(forest)
        quadrants_wrap = unsafe_wrap(Vector{Ptr{p8est_quadrant_t}}, quadrants, 2^DIM)
        levels = Vector{Int}(undef, 2^DIM)
        ps_datas = Vector{PsData}(undef, 2^DIM)
        for i = 1:(2^DIM)
            qp = PointerWrapper(quadrants_wrap[i])
            dp = PointerWrapper(P4estPsData, qp.p.user_data[])
            ps_data = unsafe_pointer_to_objref(pointer(dp.ps_data))
            isa(ps_data, InsideSolidData)&&return Cint(0)
            levels[i] = qp.level[]
            ps_datas[i] = ps_data
        end
        ka = unsafe_pointer_to_objref(pointer(fp.user_pointer))
        ps_coarsen_flag(ps_datas, levels, ka)
    end
end
"""
$(TYPEDSIGNATURES)
"""
function ps_coarsen!(p4est::Ptr{p4est_t}; recursive = 0)
    p4est_coarsen_ext(
        p4est,
        recursive,
        0,
        @cfunction(
            ps_coarsen_flag,
            Cint,
            (Ptr{p4est_t}, p4est_topidx_t, Ptr{Ptr{p4est_quadrant_t}})
        ),
        C_NULL,
        @cfunction(
            p4est_replace,
            Cvoid,
            (
                Ptr{p4est_t},
                p4est_topidx_t,
                Cint,
                Ptr{Ptr{p4est_quadrant_t}},
                Cint,
                Ptr{Ptr{p4est_quadrant_t}},
            )
        )
    )
end
"""
$(TYPEDSIGNATURES)
"""
function ps_coarsen!(p4est::Ptr{p8est_t}; recursive = 0)
    p8est_coarsen_ext(
        p4est,
        recursive,
        0,
        @cfunction(
            ps_coarsen_flag,
            Cint,
            (Ptr{p8est_t}, p4est_topidx_t, Ptr{Ptr{p8est_quadrant_t}})
        ),
        C_NULL,
        @cfunction(
            p4est_replace,
            Cvoid,
            (
                Ptr{p8est_t},
                p4est_topidx_t,
                Cint,
                Ptr{Ptr{p8est_quadrant_t}},
                Cint,
                Ptr{Ptr{p8est_quadrant_t}},
            )
        )
    )
end
"""
$(TYPEDSIGNATURES)
"""
function ps_balance!(p4est::Ptr{p4est_t})
    p4est_balance_ext(
        p4est,
        P4EST_CONNECT_FULL,
        C_NULL,
        @cfunction(
            p4est_replace,
            Cvoid,
            (
                Ptr{p4est_t},
                p4est_topidx_t,
                Cint,
                Ptr{Ptr{p4est_quadrant_t}},
                Cint,
                Ptr{Ptr{p4est_quadrant_t}},
            )
        )
    )
end
function ps_balance!(p4est::Ptr{p8est_t})
    p8est_balance_ext(
        p4est,
        P8EST_CONNECT_FULL,
        C_NULL,
        @cfunction(
            p4est_replace,
            Cvoid,
            (
                Ptr{p8est_t},
                p4est_topidx_t,
                Cint,
                Ptr{Ptr{p8est_quadrant_t}},
                Cint,
                Ptr{Ptr{p8est_quadrant_t}},
            )
        )
    )
end
function pre_ps_balance!(p4est::Ptr{p4est_t})
    p4est_balance_ext(
        p4est,
        P4EST_CONNECT_FULL,
        C_NULL,
        @cfunction(
            pre_ps_replace,
            Cvoid,
            (
                Ptr{p4est_t},
                p4est_topidx_t,
                Cint,
                Ptr{Ptr{p4est_quadrant_t}},
                Cint,
                Ptr{Ptr{p4est_quadrant_t}},
            )
        )
    )
end
function pre_ps_balance!(p4est::Ptr{p8est_t})
    p8est_balance_ext(
        p4est,
        P8EST_CONNECT_FULL,
        C_NULL,
        @cfunction(
            pre_ps_replace,
            Cvoid,
            (
                Ptr{p8est_t},
                p4est_topidx_t,
                Cint,
                Ptr{Ptr{p8est_quadrant_t}},
                Cint,
                Ptr{Ptr{p8est_quadrant_t}},
            )
        )
    )
end

"""
$(TYPEDSIGNATURES)
Update the globally maximum gradients of conserved variables as a referrence of the AMR in physical space.
"""
function update_gradmax!(ka::KA{DIM}) where {DIM}
    gradmax = ka.kinfo.status.gradmax
    gradmax .= 0.0
    trees = ka.kdata.field.trees.data
    for tree in trees
        for ps_data in tree
            isa(ps_data, InsideSolidData)&&continue
            ps_data.bound_enc<0&&continue
            for i = 1:DIM
                for j in eachindex(gradmax)
                    gradmax[j] = max(gradmax[j], abs(ps_data.sw[j, i]))
                end
            end
        end
    end
    gradmax .= MPI.Allreduce(gradmax, MPI.MAX, MPI.COMM_WORLD)
    return nothing
end

"""
$(TYPEDSIGNATURES)
Run one physical-space adaptation pass: recompute the slopes and the refinement criterion, then
refine, coarsen and 2:1-balance the forest. With `recursive = true` a single call refines /
coarsens recursively up to the level limits.

This low-level pass does **not** rebuild the ghost/neighbor/face data structures — call
[`amr_recover!`](@ref) afterwards. For routine use prefer the higher-level
[`adaptive_mesh_refinement!`](@ref) (interval-driven, recovers automatically) or simply
[`solve!`](@ref), which also uses it for the optional initial pre-refinement.
"""
function ps_adaptive_mesh_refinement!(p4est::P_pxest_t, ka::KA; recursive = false)
    # update_gradmax!(ka)
    slope!(p4est, ka)
    update_criterion!(ka)
    ps_refine!(p4est, ka; recursive = recursive ? 1 : 0)
    ps_coarsen!(p4est; recursive = recursive ? 1 : 0)
    ps_balance!(p4est)
    return nothing
end
