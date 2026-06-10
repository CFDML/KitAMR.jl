
"""
$(TYPEDSIGNATURES)
Fill the reused scratch vectors `cdf_i` (length `NDF`) and `mid_i` (length `DIM`) with the
criterion distribution `df + max_d|sdf*ds|` and the midpoint of velocity cell `c`, so the
refine/coarsen decision can evaluate the contribution criteria per cell without allocating a
full criterion matrix.  `ds` is the *physical* cell size `ps_data.ds`.
"""
@inline function _criterion_cell!(
    cdf_i,
    mid_i,
    vs_data::AbstractVsData{DIM,NDF},
    ds,
    c::Int,
) where {DIM,NDF}
    df = vs_data.df;
    sdf = vs_data.sdf
    @inbounds for d = 1:DIM
        mid_i[d] = vs_data.midpoint[c, d]
    end
    @inbounds for k = 1:NDF
        m = 0.0
        for d = 1:DIM
            a = abs(sdf[c, k, d] * ds[d]);
            a > m && (m = a)
        end
        cdf_i[k] = df[c, k] + m
    end
end

"""
$(TYPEDSIGNATURES)
"""
function vs_refine!(va_data::Velocity_Adaptive_Data, ka::KA{DIM,NDF}) where {DIM,NDF}
    trees = ka.kdata.field.trees;
    kinfo = ka.kinfo
    !isa(kinfo.config.quadrature, Vector)&&return nothing
    ds = [
        (kinfo.config.quadrature[2*i] - kinfo.config.quadrature[2*i-1]) /
        kinfo.config.vs_trees_num[i] for i = 1:DIM
    ]
    va_flags = va_data.va_flags
    id = 0
    lohner = kinfo.config.solver.ADAPT_VS_MODE === :lohner
    vmin = ntuple(d -> kinfo.config.quadrature[2*d-1], DIM)
    ds0 = ntuple(d -> ds[d], DIM)
    vstn = kinfo.config.vs_trees_num
    maxlevel = kinfo.config.solver.AMR_VS_MAXLEVEL
    τL = kinfo.config.solver.ADAPT_COEFFI_VS_LOHNER
    vsidx = lohner ? VsNeighborIndex{DIM}() : nothing
    refine_flags = Bool[]
    cdf_i = Vector{Float64}(undef, NDF)
    mid_i = Vector{Float64}(undef, DIM)
    for i in eachindex(trees.data)
        for j in eachindex(trees.data[i])
            id += 1
            ps_data = trees.data[i][j]
            isa(ps_data, InsideSolidData) && continue
            vs_data = ps_data.vs_data
            U = @view(ps_data.prim[2:(1+DIM)])
            n = vs_data.vs_num
            s1 = 0.0;
            s2 = 0.0
            if lohner
                build_vs_index!(vsidx, vs_data, vmin, ds0, vstn, maxlevel)
                s1, s2 = vs_lohner_scales(vs_data)
            end
            resize!(refine_flags, n)
            @inbounds for c = 1:n
                _criterion_cell!(cdf_i, mid_i, vs_data, ps_data.ds, c)
                base =
                    lohner ?
                    (
                        vs_lohner_indicator(vs_data, vsidx, c, s1, s2) > τL ||
                        local_contribution_refine_flag(
                            ps_data.w,
                            U,
                            mid_i,
                            cdf_i,
                            vs_data.weight[c],
                            kinfo,
                        )
                    ) :
                    contribution_refine_flag(
                        ps_data.w,
                        U,
                        mid_i,
                        cdf_i,
                        vs_data.weight[c],
                        va_data.vr,
                        kinfo,
                    )
                refine_flags[c] = vs_data.level[c] < maxlevel && base
            end
            refine_grid_stream!(vs_data, refine_flags, ds) && (va_flags[id] = true)
        end
    end
    return nothing
end
"""
$(TYPEDSIGNATURES)
"""
function vs_coarsen!(va_data::Velocity_Adaptive_Data, ka::KA{DIM,NDF}) where {DIM,NDF}
    trees = ka.kdata.field.trees;
    kinfo = ka.kinfo
    ds = [
        (kinfo.config.quadrature[2*i] - kinfo.config.quadrature[2*i-1]) /
        kinfo.config.vs_trees_num[i] for i = 1:DIM
    ]
    va_flags = va_data.va_flags
    id = 0
    lohner = kinfo.config.solver.ADAPT_VS_MODE === :lohner
    vmin = ntuple(d -> kinfo.config.quadrature[2*d-1], DIM)
    ds0 = ntuple(d -> ds[d], DIM)
    vstn = kinfo.config.vs_trees_num
    maxlevel = kinfo.config.solver.AMR_VS_MAXLEVEL
    τc = VS_LOHNER_COARSEN_RATIO * kinfo.config.solver.ADAPT_COEFFI_VS_LOHNER
    vsidx = lohner ? VsNeighborIndex{DIM}() : nothing
    coarsen_ok = Bool[]
    cdf_i = Vector{Float64}(undef, NDF)
    mid_i = Vector{Float64}(undef, DIM)
    for i in eachindex(trees.data)
        for j in eachindex(trees.data[i])
            id += 1
            ps_data = trees.data[i][j]
            isa(ps_data, InsideSolidData) && continue
            vs_data = ps_data.vs_data
            U = @view(ps_data.prim[2:(1+DIM)])
            n = vs_data.vs_num
            s1 = 0.0;
            s2 = 0.0
            if lohner
                build_vs_index!(vsidx, vs_data, vmin, ds0, vstn, maxlevel)
                s1, s2 = vs_lohner_scales(vs_data)
            end
            resize!(coarsen_ok, n)
            @inbounds for c = 1:n
                _criterion_cell!(cdf_i, mid_i, vs_data, ps_data.ds, c)
                coarsen_ok[c] =
                    lohner ?
                    (
                        vs_lohner_indicator(vs_data, vsidx, c, s1, s2) < τc &&
                        local_contribution_coarsen_flag(
                            ps_data.w,
                            U,
                            mid_i,
                            cdf_i,
                            vs_data.weight[c],
                            kinfo,
                        )
                    ) :
                    (
                        local_contribution_coarsen_flag(
                            ps_data.w,
                            U,
                            mid_i,
                            cdf_i,
                            vs_data.weight[c],
                            kinfo,
                        ) && global_contribution_coarsen_flag(
                            U,
                            mid_i,
                            cdf_i,
                            vs_data.weight[c],
                            va_data.vr,
                            kinfo,
                        )
                    )
            end
            coarsen_grid_stream!(vs_data, coarsen_ok, ds, maxlevel) && (va_flags[id] = true)
        end
    end
    return nothing
end

"""
$(TYPEDSIGNATURES)
"""
function vs_conserved_correction!(va_data::Velocity_Adaptive_Data, ka)
    trees = ka.kdata.field.trees
    va_flags = va_data.va_flags
    id = 0
    @inbounds for i in eachindex(trees.data)
        for j in eachindex(trees.data[i])
            id+=1;
            !va_flags[id]&&continue
            ps_data = trees.data[i][j]
            (isa(ps_data, InsideSolidData)||ps_data.bound_enc<0)&&continue
            vs_data = ps_data.vs_data
            conserved_I_porjection!(vs_data, ps_data.w)
        end
    end
end

function vs_resolution(ka)
    trees = ka.kdata.field.trees
    vs_resolution(trees, ka.kinfo)
end
function vs_resolution(trees::PsTrees, kinfo)
    density_res = 0.0;
    energy_res = 0.0
    @inbounds for tree in trees.data
        for ps_data in tree
            (isa(ps_data, InsideSolidData)||ps_data.bound_enc<0)&&continue
            density_res_i, energy_res_i = vs_resolution(ps_data, kinfo)
            density_res = max(density_res_i, density_res)
            energy_res = max(energy_res_i, energy_res)
        end
    end
    density_res = MPI.Allreduce(density_res, MPI.MAX, MPI.COMM_WORLD)
    energy_res = MPI.Allreduce(energy_res, MPI.MAX, MPI.COMM_WORLD)
    return Velocity_Resolution(density_res, energy_res)
end
function vs_resolution(ps_data::PsData{DIM,NDF}, kinfo) where {DIM,NDF}
    vs_data = ps_data.vs_data
    U = ps_data.prim[2:(DIM+1)]
    density_max = maximum(vs_data.df)
    c2 = @views [sum((U-u) .^ 2) for u in eachrow(vs_data.midpoint)]
    energy_max =
        NDF==1 ? 0.5*maximum(vs_data.df .* c2) :
        0.5*maximum(@view(vs_data.df[:, 1]) .* c2+@view(vs_data.df[:, 2]))
    vs_trees_num = kinfo.config.vs_trees_num
    quadrature = kinfo.config.quadrature
    du = [quadrature[2*i]-quadrature[2*i-1] for i = 1:DIM]
    weight =
        reduce(*, du)/reduce(*, vs_trees_num)/2^(DIM*(kinfo.config.solver.AMR_VS_MAXLEVEL))
    return density_max*weight, energy_max*weight
end
"""
$(TYPEDSIGNATURES)
"""
function vs_adaptive_mesh_refinement!(ka; vs_balance = false)
    vr = vs_resolution(ka)
    fp = PointerWrapper(ka.kinfo.forest.p4est)
    va_flags = zeros(Bool, fp.local_num_quadrants[])
    va_data = Velocity_Adaptive_Data(vr, va_flags)
    vs_refine!(va_data, ka)
    vs_coarsen!(va_data, ka)
    changed = any(va_flags)
    if vs_balance
        changed |= vs_balance!(ka)
    end
    vs_conserved_correction!(va_data, ka)
    return Bool(MPI.Allreduce(Int(changed), +, MPI.COMM_WORLD) > 0)
end

function initial_vs_adaptive_mesh_refinement!(
    prim,
    vs_data,
    kinfo::KInfo{DIM,NDF},
) where {DIM,NDF}
    ds = [
        (kinfo.config.quadrature[2*i] - kinfo.config.quadrature[2*i-1]) /
        kinfo.config.vs_trees_num[i] for i = 1:DIM
    ]
    maxlevel = kinfo.config.solver.AMR_VS_MAXLEVEL
    ddus = [ds ./ 2.0^i for i = 0:maxlevel]      # cell size per level
    U = @view(prim[2:(1+DIM)])
    n = vs_data.vs_num
    refine_flags = Vector{Bool}(undef, n)
    I0buf = Vector{Float64}(undef, DIM)
    I2buf = Vector{Float64}(undef, DIM)
    @inbounds for c = 1:n
        L = vs_data.level[c]
        mid = @view(vs_data.midpoint[c, :])
        refine_flags[c] =
            L < maxlevel &&
            maxwellian_refine_flag(mid, ddus[L+1], U, prim, kinfo, I0buf, I2buf)
    end
    refine_grid_stream!(vs_data, refine_flags, ds)
    return nothing
end

