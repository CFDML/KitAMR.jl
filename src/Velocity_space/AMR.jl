
"""
$(TYPEDSIGNATURES)
Fill the reused scratch vectors `cdf_i` (length `NDF`) and `mid_i` (length `DIM`) with the
criterion distribution `df + max_d|sdf*ds|` and the midpoint of velocity cell `c`, so the
refine/coarsen decision can evaluate the contribution criteria per cell without allocating a
full criterion matrix.  `ds` is the *physical* cell size `ps_data.ds`.
"""
@inline function _criterion_cell!(cdf_i, mid_i, vs_data::AbstractVsData{DIM,NDF}, ds, c::Int) where {DIM,NDF}
    df = vs_data.df; sdf = vs_data.sdf
    @inbounds for d in 1:DIM
        mid_i[d] = vs_data.midpoint[c, d]
    end
    @inbounds for k in 1:NDF
        m = 0.0
        for d in 1:DIM
            a = abs(sdf[c, k, d] * ds[d]); a > m && (m = a)
        end
        cdf_i[k] = df[c, k] + m
    end
end

const VS_HAAR_HEATFLUX_GUARD_RATIO = 4.0

@inline function _vs_effective_maxlevel(vs_data::VsData, kinfo::KInfo)
    solver = kinfo.config.solver
    solver.AMR_VS_LOCAL_LMAX || return solver.AMR_VS_MAXLEVEL
    return clamp(Int(vs_data.local_maxlevel), _vs_local_lmax_floor(solver.AMR_VS_MAXLEVEL),
                 solver.AMR_VS_MAXLEVEL)
end

function initialize_vs_local_maxlevel!(vs_data::VsData, prims, kinfo::KInfo)
    maxlevel = kinfo.config.solver.AMR_VS_MAXLEVEL
    vs_data.local_maxlevel = Int8(maxlevel)
    (kinfo.config.solver.AMR_VS_LOCAL_LMAX && isa(kinfo.config.quadrature, Vector)) ||
        return vs_data
    local_maxlevel = _vs_local_lmax_floor(maxlevel)
    for prim in prims
        local_maxlevel = max(local_maxlevel, analytic_maxwellian_local_lmax(prim, kinfo))
    end
    vs_data.local_maxlevel = Int8(clamp(local_maxlevel, _vs_local_lmax_floor(maxlevel), maxlevel))
    return vs_data
end

function adjust_local_vs_maxlevel!(ps_data::PsData{DIM,NDF}, kinfo::KInfo{DIM,NDF}, ds) where {DIM,NDF}
    vs_data = ps_data.vs_data
    solver = kinfo.config.solver
    global_maxlevel = solver.AMR_VS_MAXLEVEL
    floor_level = _vs_local_lmax_floor(global_maxlevel)
    current = clamp(Int(vs_data.local_maxlevel), floor_level, global_maxlevel)
    vs_data.local_maxlevel = Int8(current)
    global_maxlevel <= floor_level && return false

    threshold = vs_local_lmax_haar_threshold(kinfo)
    heatflux_threshold = AMR_VS_HAAR_HEATFLUX_THRESHOLD
    lower_threshold = AMR_VS_LOCAL_LMAX_COARSEN_RATIO * threshold
    lower_heatflux_threshold = AMR_VS_LOCAL_LMAX_COARSEN_RATIO * heatflux_threshold
    current_detail = vs_haar_level_max_rel_detail(vs_data, current)
    current_heatflux_detail =
        solver.AMR_VS_MODE === :haar ?
        vs_heatflux_haar_level_max_rel_detail(vs_data, ps_data.prim, current) : 0.0

    if current < global_maxlevel &&
       (current_detail > threshold || current_heatflux_detail > heatflux_threshold)
        vs_data.local_maxlevel = Int8(current + 1)
        return true
    end

    if current > floor_level &&
       current_detail < lower_threshold &&
       current_heatflux_detail < lower_heatflux_threshold
        candidate = current - 1
        coarse_detail =
            vs_haar_virtual_coarsened_level_detail(vs_data, ds, global_maxlevel, candidate)
        coarse_heatflux_detail =
            solver.AMR_VS_MODE === :haar ?
            vs_heatflux_haar_virtual_coarsened_level_detail(vs_data, ps_data.prim, ds,
                                                            global_maxlevel, candidate) : 0.0
        if coarse_detail <= threshold && coarse_heatflux_detail <= heatflux_threshold
            vs_data.local_maxlevel = Int8(candidate)
            return true
        end
    end
    return false
end

function adjust_local_vs_maxlevels!(ka::KA{DIM,NDF}) where {DIM,NDF}
    kinfo = ka.kinfo
    (kinfo.config.solver.AMR_VS_LOCAL_LMAX && isa(kinfo.config.quadrature, Vector)) ||
        return false
    ds = [(kinfo.config.quadrature[2*i] - kinfo.config.quadrature[2*i-1]) /
          kinfo.config.vs_trees_num[i] for i in 1:DIM]
    changed = false
    for tree in ka.kdata.field.trees.data
        for ps_data in tree
            isa(ps_data, InsideSolidData) && continue
            ps_data.bound_enc < 0 && continue
            adjust_local_vs_maxlevel!(ps_data, kinfo, ds) && (changed = true)
        end
    end
    return Bool(MPI.Allreduce(Int(changed), +, MPI.COMM_WORLD) > 0)
end

function _vs_group_max_contribution(ps_data::PsData{DIM,NDF}, first::Integer, nc::Integer,
                                    cdf_i, mid_i, kinfo::KInfo{DIM,NDF}) where {DIM,NDF}
    vs_data = ps_data.vs_data
    U = @view(ps_data.prim[2:1+DIM])
    max_ratio = 0.0
    @inbounds for g in 0:nc-1
        c = first + g
        _criterion_cell!(cdf_i, mid_i, vs_data, ps_data.ds, c)
        ratio = local_contribution_ratio(ps_data.w, U, mid_i, cdf_i, vs_data.weight[c], kinfo)
        ratio > max_ratio && (max_ratio = ratio)
    end
    return max_ratio
end

function _vs_haar_refine_flags!(
    refine_flags::AbstractVector{Bool},
    ps_data::PsData{DIM,NDF},
    kinfo::KInfo{DIM,NDF},
    local_maxlevel::Integer,
    detail_threshold::Real,
    heatflux_threshold::Real,
    contribution_floor::Real,
    cdf_i,
    mid_i,
) where {DIM,NDF}
    vs_data = ps_data.vs_data
    fill!(refine_flags, false)
    local_maxlevel <= 0 && return refine_flags
    nc = 2^DIM
    heatflux_contribution_floor = VS_HAAR_HEATFLUX_GUARD_RATIO * contribution_floor
    walk_maxlevel = max(Int(local_maxlevel), maximum(Int.(vs_data.level)))
    scales = _vs_component_scales(vs_data)
    U = @view(ps_data.prim[2:1+DIM])
    heatflux_scales = _vs_heatflux_scales(vs_data, U)

    @inbounds for c in 1:vs_data.vs_num
        Int(vs_data.level[c]) < local_maxlevel || continue
        _criterion_cell!(cdf_i, mid_i, vs_data, ps_data.ds, c)
        ratio =
            local_heatflux_contribution_ratio(ps_data.w, U, mid_i, cdf_i, vs_data.weight[c], kinfo)
        refine_flags[c] = ratio > heatflux_contribution_floor
    end

    for (first, level) in _vs_sibling_groups(vs_data, walk_maxlevel)
        level < local_maxlevel || continue
        df_detail = _vs_haar_group_rel_detail(vs_data, first, scales)
        heatflux_detail =
            _vs_heatflux_haar_group_rel_detail(vs_data, first, heatflux_scales, U)
        (df_detail > detail_threshold || heatflux_detail > heatflux_threshold) || continue
        @inbounds for g in 0:nc-1
            refine_flags[first + g] = true
        end
    end
    return refine_flags
end

@inline function _vs_virtual_same_level(levels::AbstractVector, first::Int, level::Int,
                                       n::Int, nc::Int)
    first + nc - 1 <= n || return false
    @inbounds for g in 1:nc-1
        Int(levels[first + g]) == level || return false
    end
    return true
end

@inline function _vs_virtual_heatflux_energy_density(vdf::AbstractMatrix, i::Integer,
                                                     c2::Real, ::Val{NDF}) where {NDF}
    if NDF == 2
        return 0.5 * (c2 * vdf[i, 1] + vdf[i, 2])
    else
        return 0.5 * c2 * vdf[i, 1]
    end
end

@inline function _vs_virtual_heatflux_integrand(vmid::AbstractMatrix, vdf::AbstractMatrix,
                                                i::Integer, U::AbstractVector, dir::Integer,
                                                ::Val{DIM}, ndf::Val{NDF}) where {DIM,NDF}
    cdir = 0.0
    c2 = 0.0
    @inbounds for d in 1:DIM
        c = vmid[i, d] - U[d]
        d == dir && (cdir = c)
        c2 += c * c
    end
    return cdir * _vs_virtual_heatflux_energy_density(vdf, i, c2, ndf)
end

function _vs_virtual_component_scales(vdf::AbstractMatrix, n::Integer,
                                      ::Val{NDF}) where {NDF}
    scales = Vector{Float64}(undef, NDF)
    @inbounds for k in 1:NDF
        s = 0.0
        for i in 1:n
            a = abs(vdf[i, k])
            a > s && (s = a)
        end
        scales[k] = max(s, EPS)
    end
    return scales
end

function _vs_virtual_heatflux_scales(vmid::AbstractMatrix, vdf::AbstractMatrix, n::Integer,
                                     U::AbstractVector, ::Val{DIM},
                                     ndf::Val{NDF}) where {DIM,NDF}
    scales = Vector{Float64}(undef, DIM)
    @inbounds for d in 1:DIM
        s = 0.0
        for i in 1:n
            a = abs(_vs_virtual_heatflux_integrand(vmid, vdf, i, U, d, Val(DIM), ndf))
            a > s && (s = a)
        end
        scales[d] = max(s, EPS)
    end
    return scales
end

function _vs_virtual_haar_group_rel_detail(vdf::AbstractMatrix, first::Int,
                                           scales::AbstractVector,
                                           ::Val{DIM}, ::Val{NDF}) where {DIM,NDF}
    nc = 2^DIM
    max_rel = 0.0
    @inbounds for k in 1:NDF
        sumsq = 0.0
        for mask in 1:nc-1
            coeff = 0.0
            for child in 1:nc
                coeff += _haar_child_sign(Val(DIM), child, mask) *
                         vdf[first + child - 1, k]
            end
            coeff /= nc
            sumsq += coeff * coeff
        end
        rel = sqrt(sumsq) / scales[k]
        rel > max_rel && (max_rel = rel)
    end
    return max_rel
end

function _vs_virtual_heatflux_haar_group_rel_detail(vmid::AbstractMatrix,
                                                    vdf::AbstractMatrix,
                                                    first::Int,
                                                    scales::AbstractVector,
                                                    U::AbstractVector,
                                                    ::Val{DIM},
                                                    ndf::Val{NDF}) where {DIM,NDF}
    nc = 2^DIM
    max_rel = 0.0
    @inbounds for d in 1:DIM
        sumsq = 0.0
        for mask in 1:nc-1
            coeff = 0.0
            for child in 1:nc
                coeff += _haar_child_sign(Val(DIM), child, mask) *
                         _vs_virtual_heatflux_integrand(vmid, vdf, first + child - 1,
                                                        U, d, Val(DIM), ndf)
            end
            coeff /= nc
            sumsq += coeff * coeff
        end
        rel = sqrt(sumsq) / scales[d]
        rel > max_rel && (max_rel = rel)
    end
    return max_rel
end

@inline function _vs_cancel_virtual_coarsen!(coarsen_ok::AbstractVector{Bool},
                                             origin_first::AbstractVector{Int},
                                             origin_count::AbstractVector{Int},
                                             i::Integer)
    n = origin_count[i]
    n > 1 || return nothing
    first = origin_first[i]
    @inbounds for g in 0:n-1
        coarsen_ok[first + g] = false
    end
    return nothing
end

function _vs_haar_coarsen_consistency_guard!(
    coarsen_ok::AbstractVector{Bool},
    ps_data::PsData{DIM,NDF},
    kinfo::KInfo{DIM,NDF},
    local_maxlevel::Integer,
    maxlevel::Integer,
    detail_threshold::Real,
    heatflux_threshold::Real,
    contribution_floor::Real,
    vs_ds,
    cdf_i,
    mid_i,
) where {DIM,NDF}
    local_maxlevel <= 0 && return coarsen_ok
    vs_data = ps_data.vs_data
    n = vs_data.vs_num
    nc = 2^DIM
    vlevel = Vector{Int8}(undef, n)
    vweight = Vector{Float64}(undef, n)
    vmid = Matrix{Float64}(undef, n, DIM)
    vdf = Matrix{Float64}(undef, n, NDF)
    origin_first = Vector{Int}(undef, n)
    origin_count = Vector{Int}(undef, n)
    merged_candidate = Vector{Bool}(undef, n)
    flag = zeros(Float64, maxlevel)
    any_candidate = false
    w = 0
    index = 1

    @inbounds while index <= n
        fl = Int(vs_data.level[index])
        if fl > 0
            aligned = flag[fl] % 1 == 0.0
            if aligned && index + nc - 1 <= n && _all_same_level(vs_data, index, vs_data.level[index], nc)
                grp = true
                for g in 0:nc-1
                    coarsen_ok[index + g] || (grp = false; break)
                end
                if grp
                    w += 1
                    parent_level = fl - 1
                    vlevel[w] = Int8(parent_level)
                    vweight[w] = vs_data.weight[index] * nc
                    for d in 1:DIM
                        vmid[w, d] =
                            vs_data.midpoint[index, d] -
                            0.5 * (vs_ds[d] / 2^fl) * RMT[DIM][1][d]
                    end
                    for k in 1:NDF
                        s = 0.0
                        for g in 0:nc-1
                            s += vs_data.df[index + g, k]
                        end
                        vdf[w, k] = s / nc
                    end
                    origin_first[w] = index
                    origin_count[w] = nc
                    merged_candidate[w] = parent_level < local_maxlevel
                    any_candidate |= merged_candidate[w]
                else
                    for g in 0:nc-1
                        w += 1
                        src = index + g
                        vlevel[w] = vs_data.level[src]
                        vweight[w] = vs_data.weight[src]
                        for d in 1:DIM
                            vmid[w, d] = vs_data.midpoint[src, d]
                        end
                        for k in 1:NDF
                            vdf[w, k] = vs_data.df[src, k]
                        end
                        origin_first[w] = src
                        origin_count[w] = 1
                        merged_candidate[w] = false
                    end
                end
                index += nc
                if fl > 1
                    for l in 1:fl-1
                        flag[l] += 1 / 2^(DIM * (fl - l))
                    end
                end
            else
                w += 1
                vlevel[w] = vs_data.level[index]
                vweight[w] = vs_data.weight[index]
                for d in 1:DIM
                    vmid[w, d] = vs_data.midpoint[index, d]
                end
                for k in 1:NDF
                    vdf[w, k] = vs_data.df[index, k]
                end
                origin_first[w] = index
                origin_count[w] = 1
                merged_candidate[w] = false
                for l in 1:fl
                    flag[l] += 1 / 2^(DIM * (fl - l + 1))
                end
                index += 1
            end
        else
            w += 1
            vlevel[w] = vs_data.level[index]
            vweight[w] = vs_data.weight[index]
            for d in 1:DIM
                vmid[w, d] = vs_data.midpoint[index, d]
            end
            for k in 1:NDF
                vdf[w, k] = vs_data.df[index, k]
            end
            origin_first[w] = index
            origin_count[w] = 1
            merged_candidate[w] = false
            index += 1
        end
    end
    any_candidate || return coarsen_ok

    U = @view(ps_data.prim[2:1+DIM])
    heatflux_contribution_floor = VS_HAAR_HEATFLUX_GUARD_RATIO * contribution_floor
    @inbounds for i in 1:w
        merged_candidate[i] || continue
        for d in 1:DIM
            mid_i[d] = vmid[i, d]
        end
        for k in 1:NDF
            cdf_i[k] = vdf[i, k]
        end
        ratio =
            local_heatflux_contribution_ratio(ps_data.w, U, mid_i, cdf_i, vweight[i], kinfo)
        ratio > heatflux_contribution_floor &&
            _vs_cancel_virtual_coarsen!(coarsen_ok, origin_first, origin_count, i)
    end

    scales = _vs_virtual_component_scales(vdf, w, Val(NDF))
    heatflux_scales = _vs_virtual_heatflux_scales(vmid, vdf, w, U, Val(DIM), Val(NDF))
    fill!(flag, 0.0)
    index = 1
    @inbounds while index <= w
        level = Int(vlevel[index])
        if level > 0
            aligned = abs(mod(flag[level], 1.0)) < 1.0e-12
            if aligned && _vs_virtual_same_level(vlevel, index, level, w, nc)
                if level < local_maxlevel
                    has_candidate = false
                    for g in 0:nc-1
                        merged_candidate[index + g] && (has_candidate = true; break)
                    end
                    if has_candidate
                        rel = _vs_virtual_haar_group_rel_detail(vdf, index, scales,
                                                                Val(DIM), Val(NDF))
                        heatflux_rel =
                            _vs_virtual_heatflux_haar_group_rel_detail(vmid, vdf, index,
                                                                       heatflux_scales, U,
                                                                       Val(DIM), Val(NDF))
                        if rel > detail_threshold || heatflux_rel > heatflux_threshold
                            for g in 0:nc-1
                                merged_candidate[index + g] &&
                                    _vs_cancel_virtual_coarsen!(coarsen_ok, origin_first,
                                                                origin_count, index + g)
                            end
                        end
                    end
                end
                index += nc
                if level > 1
                    for l in 1:level-1
                        flag[l] += 1 / 2^(DIM * (level - l))
                    end
                end
            else
                for l in 1:level
                    flag[l] += 1 / 2^(DIM * (level - l + 1))
                end
                index += 1
            end
        else
            index += 1
        end
    end
    return coarsen_ok
end

function _vs_haar_coarsen_flags!(
    coarsen_ok::AbstractVector{Bool},
    ps_data::PsData{DIM,NDF},
    kinfo::KInfo{DIM,NDF},
    local_maxlevel::Integer,
    maxlevel::Integer,
    detail_threshold::Real,
    heatflux_threshold::Real,
    contribution_floor::Real,
    vs_ds,
    cdf_i,
    mid_i,
) where {DIM,NDF}
    vs_data = ps_data.vs_data
    fill!(coarsen_ok, false)
    nc = 2^DIM
    lower_detail_threshold = AMR_VS_LOCAL_LMAX_COARSEN_RATIO * detail_threshold
    lower_heatflux_threshold = AMR_VS_LOCAL_LMAX_COARSEN_RATIO * heatflux_threshold
    coarsen_contribution_floor =
        VS_LSR_COARSEN_RATIO * VS_HAAR_HEATFLUX_GUARD_RATIO * contribution_floor
    U = @view(ps_data.prim[2:1+DIM])
    @inbounds for c in 1:vs_data.vs_num
        if Int(vs_data.level[c]) > local_maxlevel
            coarsen_ok[c] = true
            continue
        end
        _criterion_cell!(cdf_i, mid_i, vs_data, ps_data.ds, c)
        ratio =
            local_heatflux_contribution_ratio(ps_data.w, U, mid_i, cdf_i, vs_data.weight[c], kinfo)
        coarsen_ok[c] = ratio < coarsen_contribution_floor
    end
    maxlevel <= 0 && return coarsen_ok
    walk_maxlevel = max(Int(maxlevel), maximum(Int.(vs_data.level)))
    scales = _vs_component_scales(vs_data)
    heatflux_scales = _vs_heatflux_scales(vs_data, U)
    for (first, level) in _vs_sibling_groups(vs_data, walk_maxlevel)
        0 < level <= local_maxlevel || continue
        rel = _vs_haar_group_rel_detail(vs_data, first, scales)
        heatflux_rel =
            _vs_heatflux_haar_group_rel_detail(vs_data, first, heatflux_scales, U)
        (rel < lower_detail_threshold && heatflux_rel < lower_heatflux_threshold) && continue
        @inbounds for g in 0:nc-1
            coarsen_ok[first + g] = false
        end
    end
    _vs_haar_coarsen_consistency_guard!(coarsen_ok, ps_data, kinfo, local_maxlevel,
                                        maxlevel, detail_threshold, heatflux_threshold,
                                        contribution_floor, vs_ds, cdf_i, mid_i)
    return coarsen_ok
end

"""
$(TYPEDSIGNATURES)
"""
function vs_refine!(
    va_data::Velocity_Adaptive_Data,
    ka::KA{DIM,NDF};
    linear_reconstruction::Bool = false,
) where{DIM,NDF}
    trees = ka.kdata.field.trees;kinfo = ka.kinfo
    !isa(kinfo.config.quadrature,Vector)&&return nothing
    ds = [(kinfo.config.quadrature[2*i] - kinfo.config.quadrature[2*i-1]) /
    kinfo.config.vs_trees_num[i] for i in 1:DIM]
    va_flags = va_data.va_flags
    id = 0
    vs_mode = kinfo.config.solver.AMR_VS_MODE
    lsr_mode = vs_mode === :lsr
    haar_mode = vs_mode === :haar
    use_linear_reconstruction = linear_reconstruction && lsr_mode
    if lsr_mode
        vmin = ntuple(d -> kinfo.config.quadrature[2*d-1], DIM)
        ds0 = ntuple(d -> ds[d], DIM)
        vstn = kinfo.config.vs_trees_num
        vsidx = VsNeighborIndex{DIM}()
    end
    maxlevel = kinfo.config.solver.AMR_VS_MAXLEVEL
    τR = kinfo.config.solver.AMR_VS_LSR_THRES
    τH = haar_mode ? vs_local_lmax_haar_threshold(kinfo) : 0.0
    τQ = haar_mode ? AMR_VS_HAAR_HEATFLUX_THRESHOLD : 0.0
    contribution_threshold = kinfo.config.solver.AMR_VS_CONTRI_FLOOR
    refine_flags = Bool[]
    cdf_i = Vector{Float64}(undef, NDF)
    mid_i = Vector{Float64}(undef, DIM)
    if lsr_mode
        lsr_neighbors = Int[]
        lsr_normal = zeros(Float64, DIM, DIM)
        lsr_rhs = zeros(Float64, DIM)
        lsr_x = zeros(Float64, DIM)
    end
    for i in eachindex(trees.data)
        for j in eachindex(trees.data[i])
            id += 1
            ps_data = trees.data[i][j]
            isa(ps_data,InsideSolidData) && continue
            vs_data = ps_data.vs_data
            local_maxlevel = _vs_effective_maxlevel(vs_data, kinfo)
            U = @view(ps_data.prim[2:1+DIM])
            n = vs_data.vs_num
            if lsr_mode
                build_vs_index!(vsidx, vs_data, vmin, ds0, vstn, maxlevel)
                s1, s2 = vs_lsr_scales(vs_data)
            end
            resize!(refine_flags, n)
            if haar_mode
                _vs_haar_refine_flags!(refine_flags, ps_data, kinfo, local_maxlevel, τH, τQ,
                                       contribution_threshold, cdf_i, mid_i)
            else
                @inbounds for c in 1:n
                    if vs_data.level[c] >= local_maxlevel
                        refine_flags[c] = false
                        continue
                    end
                    _criterion_cell!(cdf_i, mid_i, vs_data, ps_data.ds, c)
                    ratio = local_contribution_ratio(ps_data.w, U, mid_i, cdf_i, vs_data.weight[c], kinfo)
                    if lsr_mode
                        refine_flags[c] =
                            ratio > contribution_threshold &&
                            vs_lsr_indicator(vs_data, vsidx, c, s1, s2, lsr_neighbors,
                                             lsr_normal, lsr_rhs, lsr_x) > τR
                    else
                        refine_flags[c] = ratio > contribution_threshold
                    end
                end
            end
            changed = use_linear_reconstruction ?
                      refine_grid_stream_linear!(vs_data, refine_flags, ds, vsidx) :
                      refine_grid_stream!(vs_data, refine_flags, ds)
            changed && (va_flags[id] = true)
        end
    end
    return nothing
end
"""
$(TYPEDSIGNATURES)
"""
function vs_coarsen!(va_data::Velocity_Adaptive_Data,ka::KA{DIM,NDF})where{DIM,NDF}
    trees = ka.kdata.field.trees;kinfo = ka.kinfo
    ds = [(kinfo.config.quadrature[2*i] - kinfo.config.quadrature[2*i-1]) /
    kinfo.config.vs_trees_num[i] for i in 1:DIM]
    va_flags = va_data.va_flags
    id = 0
    vs_mode = kinfo.config.solver.AMR_VS_MODE
    lsr_mode = vs_mode === :lsr
    haar_mode = vs_mode === :haar
    if lsr_mode
        vmin = ntuple(d -> kinfo.config.quadrature[2*d-1], DIM)
        ds0 = ntuple(d -> ds[d], DIM)
        vstn = kinfo.config.vs_trees_num
        vsidx = VsNeighborIndex{DIM}()
    end
    maxlevel = kinfo.config.solver.AMR_VS_MAXLEVEL
    τr_c = VS_LSR_COARSEN_RATIO * kinfo.config.solver.AMR_VS_LSR_THRES
    τH = haar_mode ? vs_local_lmax_haar_threshold(kinfo) : 0.0
    τQ = haar_mode ? AMR_VS_HAAR_HEATFLUX_THRESHOLD : 0.0
    contribution_coarsen_threshold =
        VS_LSR_COARSEN_RATIO * kinfo.config.solver.AMR_VS_CONTRI_FLOOR
    coarsen_ok = Bool[]
    cdf_i = Vector{Float64}(undef, NDF)
    mid_i = Vector{Float64}(undef, DIM)
    if lsr_mode
        lsr_neighbors = Int[]
        lsr_normal = zeros(Float64, DIM, DIM)
        lsr_rhs = zeros(Float64, DIM)
        lsr_x = zeros(Float64, DIM)
    end
    for i in eachindex(trees.data)
        for j in eachindex(trees.data[i])
            id += 1
            ps_data = trees.data[i][j]
            isa(ps_data,InsideSolidData) && continue
            vs_data = ps_data.vs_data
            local_maxlevel = _vs_effective_maxlevel(vs_data, kinfo)
            U = @view(ps_data.prim[2:1+DIM])
            n = vs_data.vs_num
            if lsr_mode
                build_vs_index!(vsidx, vs_data, vmin, ds0, vstn, maxlevel)
                s1, s2 = vs_lsr_scales(vs_data)
            end
            resize!(coarsen_ok, n)
            if haar_mode
                _vs_haar_coarsen_flags!(coarsen_ok, ps_data, kinfo, local_maxlevel, maxlevel, τH, τQ,
                                        kinfo.config.solver.AMR_VS_CONTRI_FLOOR, ds, cdf_i, mid_i)
            else
                @inbounds for c in 1:n
                    if Int(vs_data.level[c]) > local_maxlevel
                        coarsen_ok[c] = true
                        continue
                    end
                    _criterion_cell!(cdf_i, mid_i, vs_data, ps_data.ds, c)
                    ratio = local_contribution_ratio(ps_data.w, U, mid_i, cdf_i, vs_data.weight[c], kinfo)
                    if lsr_mode
                        if ratio >= contribution_coarsen_threshold
                            coarsen_ok[c] = false
                        else
                            coarsen_ok[c] =
                                vs_lsr_indicator(vs_data, vsidx, c, s1, s2, lsr_neighbors,
                                                 lsr_normal, lsr_rhs, lsr_x) < τr_c
                        end
                    else
                        coarsen_ok[c] = ratio < contribution_coarsen_threshold
                    end
                end
            end
            coarsen_grid_stream!(vs_data, coarsen_ok, ds, maxlevel) && (va_flags[id] = true)
        end
    end
    return nothing
end

"""
$(TYPEDSIGNATURES)
"""
function vs_conserved_correction!(va_data::Velocity_Adaptive_Data,ka)
    trees = ka.kdata.field.trees
    va_flags = va_data.va_flags
    id = 0
    @inbounds for i in eachindex(trees.data)
        for j in eachindex(trees.data[i])
            id+=1;!va_flags[id]&&continue
            ps_data = trees.data[i][j]
            (isa(ps_data,InsideSolidData)||ps_data.bound_enc<0)&&continue
            vs_data = ps_data.vs_data
            conserved_I_projection!(vs_data,ps_data.w)
        end
    end
end

"""
$(TYPEDSIGNATURES)
"""
function _vs_adaptive_mesh_refinement_result!(ka;vs_balance = false)
    fp = PointerWrapper(ka.kinfo.forest.p4est)
    va_flags = zeros(Bool,fp.local_num_quadrants[])
    va_data = Velocity_Adaptive_Data(va_flags)
    vs_refine!(va_data,ka)
    vs_coarsen!(va_data,ka)
    changed = any(va_flags)
    if vs_balance
        changed |= vs_balance!(ka)
    end
    vs_conserved_correction!(va_data,ka)
    return Bool(MPI.Allreduce(Int(changed), +, MPI.COMM_WORLD) > 0)
end

function vs_adaptive_mesh_refinement!(ka;vs_balance = false)
    return _vs_adaptive_mesh_refinement_result!(ka; vs_balance = vs_balance)
end

function initial_vs_adaptive_mesh_refinement!(prim::AbstractVector{<:Real},vs_data,kinfo::KInfo{DIM,NDF}) where{DIM,NDF}
    return initial_vs_adaptive_mesh_refinement!((prim,),vs_data,kinfo)
end
function initial_vs_adaptive_mesh_refinement!(prims,vs_data,kinfo::KInfo{DIM,NDF}) where{DIM,NDF}
    ds = [(kinfo.config.quadrature[2*i] - kinfo.config.quadrature[2*i-1]) /
        kinfo.config.vs_trees_num[i] for i in 1:DIM]
    maxlevel = kinfo.config.solver.AMR_VS_MAXLEVEL
    local_maxlevel = _vs_effective_maxlevel(vs_data, kinfo)
    ddus = [ds ./ 2.0^i for i in 0:maxlevel]      # cell size per level
    n = vs_data.vs_num
    refine_flags = Vector{Bool}(undef, n)
    I0buf = Vector{Float64}(undef, DIM)
    I2buf = Vector{Float64}(undef, DIM)
    @inbounds for c in 1:n
        L = vs_data.level[c]
        mid = @view(vs_data.midpoint[c, :])
        flag = false
        if L < local_maxlevel
            for prim in prims
                U = @view(prim[2:1+DIM])
                if maxwellian_refine_flag(mid, ddus[L+1], U, prim, kinfo, I0buf, I2buf)
                    flag = true
                    break
                end
            end
        end
        refine_flags[c] = flag
    end
    refine_grid_stream!(vs_data, refine_flags, ds)
    return nothing
end
