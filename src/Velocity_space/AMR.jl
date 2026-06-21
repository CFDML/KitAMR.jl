
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
    lsr_mode = kinfo.config.solver.AMR_VS_MODE === :lsr
    use_linear_reconstruction = linear_reconstruction && lsr_mode
    if lsr_mode
        vmin = ntuple(d -> kinfo.config.quadrature[2*d-1], DIM)
        ds0 = ntuple(d -> ds[d], DIM)
        vstn = kinfo.config.vs_trees_num
        vsidx = VsNeighborIndex{DIM}()
    end
    maxlevel = kinfo.config.solver.AMR_VS_MAXLEVEL
    τR = kinfo.config.solver.AMR_VS_LSR_THRES
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
            U = @view(ps_data.prim[2:1+DIM])
            n = vs_data.vs_num
            if lsr_mode
                build_vs_index!(vsidx, vs_data, vmin, ds0, vstn, maxlevel)
                s1, s2 = vs_lsr_scales(vs_data)
            end
            resize!(refine_flags, n)
            @inbounds for c in 1:n
                if vs_data.level[c] >= maxlevel
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
    lsr_mode = kinfo.config.solver.AMR_VS_MODE === :lsr
    if lsr_mode
        vmin = ntuple(d -> kinfo.config.quadrature[2*d-1], DIM)
        ds0 = ntuple(d -> ds[d], DIM)
        vstn = kinfo.config.vs_trees_num
        vsidx = VsNeighborIndex{DIM}()
    end
    maxlevel = kinfo.config.solver.AMR_VS_MAXLEVEL
    τr_c = VS_LSR_COARSEN_RATIO * kinfo.config.solver.AMR_VS_LSR_THRES
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
            U = @view(ps_data.prim[2:1+DIM])
            n = vs_data.vs_num
            if lsr_mode
                build_vs_index!(vsidx, vs_data, vmin, ds0, vstn, maxlevel)
                s1, s2 = vs_lsr_scales(vs_data)
            end
            resize!(coarsen_ok, n)
            @inbounds for c in 1:n
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
