
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
    vmin = ntuple(d -> kinfo.config.quadrature[2*d-1], DIM)
    ds0 = ntuple(d -> ds[d], DIM)
    vstn = kinfo.config.vs_trees_num
    maxlevel = kinfo.config.solver.AMR_VS_MAXLEVEL
    τR = kinfo.config.solver.ADAPT_COEFFI_VS_LSR
    lsr_floor = kinfo.config.solver.ADAPT_COEFFI_VS_LSR_FLOOR
    vsidx = VsNeighborIndex{DIM}()
    refine_flags = Bool[]
    cdf_i = Vector{Float64}(undef, NDF)
    mid_i = Vector{Float64}(undef, DIM)
    lsr_neighbors = Int[]
    lsr_normal = zeros(Float64, DIM, DIM)
    lsr_rhs = zeros(Float64, DIM)
    lsr_x = zeros(Float64, DIM)
    for i in eachindex(trees.data)
        for j in eachindex(trees.data[i])
            id += 1
            ps_data = trees.data[i][j]
            isa(ps_data,InsideSolidData) && continue
            vs_data = ps_data.vs_data
            U = @view(ps_data.prim[2:1+DIM])
            n = vs_data.vs_num
            s1 = 0.0; s2 = 0.0
            build_vs_index!(vsidx, vs_data, vmin, ds0, vstn, maxlevel)
            s1, s2 = vs_lsr_scales(vs_data)
            resize!(refine_flags, n)
            @inbounds for c in 1:n
                _criterion_cell!(cdf_i, mid_i, vs_data, ps_data.ds, c)
                η = vs_lsr_indicator(vs_data, vsidx, c, s1, s2, lsr_neighbors, lsr_normal, lsr_rhs, lsr_x)
                ratio = local_contribution_ratio(ps_data.w, U, mid_i, cdf_i, vs_data.weight[c], kinfo)
                base = η > τR && ratio > lsr_floor
                refine_flags[c] = vs_data.level[c] < maxlevel && base
            end
            changed = linear_reconstruction ?
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
    vmin = ntuple(d -> kinfo.config.quadrature[2*d-1], DIM)
    ds0 = ntuple(d -> ds[d], DIM)
    vstn = kinfo.config.vs_trees_num
    maxlevel = kinfo.config.solver.AMR_VS_MAXLEVEL
    τr_c = VS_LSR_COARSEN_RATIO * kinfo.config.solver.ADAPT_COEFFI_VS_LSR
    lsr_floor = kinfo.config.solver.ADAPT_COEFFI_VS_LSR_FLOOR
    vsidx = VsNeighborIndex{DIM}()
    coarsen_ok = Bool[]
    cdf_i = Vector{Float64}(undef, NDF)
    mid_i = Vector{Float64}(undef, DIM)
    lsr_neighbors = Int[]
    lsr_normal = zeros(Float64, DIM, DIM)
    lsr_rhs = zeros(Float64, DIM)
    lsr_x = zeros(Float64, DIM)
    for i in eachindex(trees.data)
        for j in eachindex(trees.data[i])
            id += 1
            ps_data = trees.data[i][j]
            isa(ps_data,InsideSolidData) && continue
            vs_data = ps_data.vs_data
            U = @view(ps_data.prim[2:1+DIM])
            n = vs_data.vs_num
            s1 = 0.0; s2 = 0.0
            build_vs_index!(vsidx, vs_data, vmin, ds0, vstn, maxlevel)
            s1, s2 = vs_lsr_scales(vs_data)
            resize!(coarsen_ok, n)
            @inbounds for c in 1:n
                _criterion_cell!(cdf_i, mid_i, vs_data, ps_data.ds, c)
                η = vs_lsr_indicator(vs_data, vsidx, c, s1, s2, lsr_neighbors, lsr_normal, lsr_rhs, lsr_x)
                ratio = local_contribution_ratio(ps_data.w, U, mid_i, cdf_i, vs_data.weight[c], kinfo)
                coarsen_ok[c] = (η < τr_c && ratio < lsr_floor) || ratio < 0.5 * lsr_floor
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
function _record_vs_counts!(counts::AbstractVector{Int}, trees::PsTrees)
    id = 0
    @inbounds for tree in trees.data
        for ps_data in tree
            id += 1
            counts[id] = isa(ps_data, InsideSolidData) ? 0 : ps_data.vs_data.vs_num
        end
    end
    return counts
end

function _max_vs_count_change_ratio(before::AbstractVector{Int}, trees::PsTrees)
    id = 0
    ratio = 0.0
    @inbounds for tree in trees.data
        for ps_data in tree
            id += 1
            isa(ps_data, InsideSolidData) && continue
            n0 = before[id]
            n0 > 0 || continue
            n1 = ps_data.vs_data.vs_num
            ratio = max(ratio, abs(n1 - n0) / n0)
        end
    end
    return ratio
end

function _vs_adaptive_mesh_refinement_result!(ka;vs_balance = false)
    trees = ka.kdata.field.trees
    fp = PointerWrapper(ka.kinfo.forest.p4est)
    before = Vector{Int}(undef, fp.local_num_quadrants[])
    _record_vs_counts!(before, trees)
    va_flags = zeros(Bool,fp.local_num_quadrants[])
    va_data = Velocity_Adaptive_Data(va_flags)
    vs_refine!(va_data,ka)
    vs_coarsen!(va_data,ka)
    changed = any(va_flags)
    if vs_balance
        changed |= vs_balance!(ka)
    end
    vs_conserved_correction!(va_data,ka)
    change_ratio = _max_vs_count_change_ratio(before, trees)
    return Bool(MPI.Allreduce(Int(changed), +, MPI.COMM_WORLD) > 0),
           MPI.Allreduce(change_ratio, MPI.MAX, MPI.COMM_WORLD)
end

function vs_adaptive_mesh_refinement!(ka;vs_balance = false)
    changed, _ = _vs_adaptive_mesh_refinement_result!(ka; vs_balance = vs_balance)
    return changed
end

function initial_vs_adaptive_mesh_refinement!(prim::AbstractVector{<:Real},vs_data,kinfo::KInfo{DIM,NDF}) where{DIM,NDF}
    return initial_vs_adaptive_mesh_refinement!((prim,),vs_data,kinfo)
end
function initial_vs_adaptive_mesh_refinement!(prims,vs_data,kinfo::KInfo{DIM,NDF}) where{DIM,NDF}
    ds = [(kinfo.config.quadrature[2*i] - kinfo.config.quadrature[2*i-1]) /
        kinfo.config.vs_trees_num[i] for i in 1:DIM]
    maxlevel = kinfo.config.solver.AMR_VS_MAXLEVEL
    ddus = [ds ./ 2.0^i for i in 0:maxlevel]      # cell size per level
    n = vs_data.vs_num
    refine_flags = Vector{Bool}(undef, n)
    I0buf = Vector{Float64}(undef, DIM)
    I2buf = Vector{Float64}(undef, DIM)
    @inbounds for c in 1:n
        L = vs_data.level[c]
        mid = @view(vs_data.midpoint[c, :])
        flag = false
        if L < maxlevel
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
