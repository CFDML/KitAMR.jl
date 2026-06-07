function pre_vs_refine!(trees::PsTrees{DIM,NDF}, kinfo::KInfo{DIM,NDF}) where{DIM,NDF}
    vr = vs_resolution(trees,kinfo)
    !isa(kinfo.config.quadrature,Vector)&&return nothing
    ds = [(kinfo.config.quadrature[2*i] - kinfo.config.quadrature[2*i-1]) /
        kinfo.config.vs_trees_num[i] for i in 1:DIM]
    vs_refine_udf = kinfo.config.user_defined.static_vs_refine_flag
    for _ = 1:kinfo.config.solver.AMR_VS_MAXLEVEL
        for i in eachindex(trees.data)
            for j in eachindex(trees.data[i])
                ps_data = trees.data[i][j]
                isa(ps_data,InsideSolidData) && continue
                vs_data = ps_data.vs_data
                U = ps_data.prim[2:1+DIM]
                lnmidpoint = reshape(vs_data.midpoint, :)
                lndf = reshape(vs_data.df, :)
                lnsdf = reshape(vs_data.sdf, :)
                lnflux = reshape(vs_data.flux, :)
                index = 1
                midpoint_index = Vector{Int}(undef, DIM)
                df_index = Vector{Int}(undef, NDF)
                while index < vs_data.vs_num + 1
                    @inbounds for i = 1:DIM
                        midpoint_index[i] = (i - 1) * (vs_data.vs_num) + index
                    end
                    @inbounds for i = 1:NDF
                        df_index[i] = (i - 1) * (vs_data.vs_num) + index
                    end
                    midpoint = @view(lnmidpoint[midpoint_index])
                    df = @view(lndf[df_index])
                    if vs_data.level[index] < kinfo.config.solver.AMR_VS_MAXLEVEL &&
                        (   vs_refine_udf==null_udf ? 
                            contribution_refine_flag(ps_data.w, U, midpoint, df, vs_data.weight[index], vr, kinfo) : 
                            vs_refine_udf(midpoint;level = vs_data.level[index],du = ds./2^vs_data.level[index])
                        )
                        midpoint_new = midpoint_refine(DIM,midpoint, vs_data.level[index], ds)
                        df_new = df_refine(DIM,midpoint, midpoint_new, df)
                        vs_data.vs_num += 2^DIM - 1
                        level_refine_replace!(DIM,vs_data.level, index)
                        weight_refine_replace!(DIM,vs_data.weight, index)
                        midpoint_refine_replace!(
                            DIM,
                            lnmidpoint,
                            midpoint_new,
                            vs_data.vs_num,
                            index,
                        )
                        df_refine_replace!(DIM,NDF,lndf, df_new, vs_data.vs_num, index)
                        sdf_refine_replace!(DIM,NDF,lnsdf)
                        flux_refine_replace!(DIM,NDF,lnflux)
                        index += 2^DIM - 1
                    end
                    index += 1
                end
                vs_data.midpoint = reshape(lnmidpoint, vs_data.vs_num, DIM)
                vs_data.df = reshape(lndf, vs_data.vs_num, NDF)
                vs_data.sdf = reshape(lnsdf, vs_data.vs_num, NDF, DIM)
                vs_data.flux = reshape(lnflux, vs_data.vs_num, NDF)
            end
        end
    end
    return nothing
end
function criterion_df(ps_data::AbstractPsData{DIM,NDF}) where{DIM,NDF}
    vs_data = ps_data.vs_data;df = vs_data.df;sdf = vs_data.sdf
    ds = ps_data.ds
    cdf = similar(vs_data.df)
    ddfi = zeros(DIM)
    for j in axes(cdf,2)
        for i in axes(cdf,1)
            for k in 1:DIM
                ddfi[k] = abs(sdf[i,j,k]*ds[k])
            end
            cdf[i,j] = df[i,j]+maximum(ddfi)
        end
    end
    return cdf
end

"""
$(TYPEDSIGNATURES)
Fill the reused scratch vectors `cdf_i` (length `NDF`) and `mid_i` (length `DIM`) with the
criterion distribution `df + max_d|sdf*ds|` and the midpoint of velocity cell `c`.  Lets the
refine/coarsen decision evaluate the contribution criteria per cell without allocating a full
`criterion_df` matrix.  `ds` is the *physical* cell size `ps_data.ds`.
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
function vs_refine!(va_data::Velocity_Adaptive_Data, ka::KA{DIM,NDF}) where{DIM,NDF}
    trees = ka.kdata.field.trees;kinfo = ka.kinfo
    !isa(kinfo.config.quadrature,Vector)&&return nothing
    ds = [(kinfo.config.quadrature[2*i] - kinfo.config.quadrature[2*i-1]) /
    kinfo.config.vs_trees_num[i] for i in 1:DIM]
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
            isa(ps_data,InsideSolidData) && continue
            vs_data = ps_data.vs_data
            U = @view(ps_data.prim[2:1+DIM])
            n = vs_data.vs_num
            s1 = 0.0; s2 = 0.0
            if lohner
                build_vs_index!(vsidx, vs_data, vmin, ds0, vstn, maxlevel)
                s1, s2 = vs_lohner_scales(vs_data)
            end
            resize!(refine_flags, n)
            @inbounds for c in 1:n
                _criterion_cell!(cdf_i, mid_i, vs_data, ps_data.ds, c)
                base = lohner ?
                    (vs_lohner_indicator(vs_data, vsidx, c, s1, s2) > τL ||
                     local_contribution_refine_flag(ps_data.w, U, mid_i, cdf_i, vs_data.weight[c], kinfo)) :
                    contribution_refine_flag(ps_data.w, U, mid_i, cdf_i, vs_data.weight[c], va_data.vr, kinfo)
                refine_flags[c] = vs_data.level[c] < maxlevel && base
            end
            refine_grid_stream!(vs_data, refine_flags, ds) && (va_flags[id] = true)
        end
    end
    return nothing
end
"""
$(SIGNATURES)
"""
function midpoint_refine_replace!(DIM::Integer,lnmidpoint::AbstractVector, midpoint_new::AbstractMatrix, vs_num::Int, index::Int)
    for i = 1:DIM
        deleteat!(lnmidpoint, (i - 1) * (vs_num) + index)
        for j = 2^DIM:-1:1
            insert!(lnmidpoint, (i - 1) * (vs_num) + index, midpoint_new[j, i])
        end
    end
end
"""
$(SIGNATURES)
"""
function weight_refine_replace!(DIM::Integer,weight::AbstractVector, index::Int)
    weight_new = popat!(weight, index) / 2^DIM
    for _ = 1:2^DIM
        insert!(weight, index, weight_new)
    end
end
"""
$(SIGNATURES)
"""
function level_refine_replace!(DIM::Integer,level::AbstractVector, index::Int)
    level_new = popat!(level, index) + 1
    for _ = 1:2^DIM
        insert!(level, index, level_new)
    end
end
"""
$(SIGNATURES)
"""
function df_refine_replace!(DIM::Integer,NDF::Integer,lndf::AbstractVector, df_new::AbstractMatrix, vs_num::Int, index::Int)
    for i = 1:NDF
        deleteat!(lndf, (i - 1) * (vs_num) + index)
        for j = 2^DIM:-1:1
            insert!(lndf, (i - 1) * (vs_num) + index, df_new[j, i])
        end
    end
end
"""
$(SIGNATURES)
"""
function sdf_refine_replace!(DIM::Integer,NDF::Integer,lnsdf::AbstractVector)
    append!(lnsdf, zeros((2^DIM - 1) * NDF * DIM))
end
"""
$(SIGNATURES)
"""
function flux_refine_replace!(DIM::Integer,NDF::Integer,lnflux::AbstractVector)
    append!(lnflux, zeros((2^DIM - 1) * NDF))
end
"""
$(SIGNATURES)
"""
function midpoint_refine(DIM::Integer,midpoint::AbstractVector, level::Int8, ds::AbstractVector)
    midpoint_new = Matrix{Float64}(undef, 2^DIM, DIM)
    ds_new = ds / 2^(level + 1)
    for i = 1:2^DIM
        midpoint_new[i, :] .= @. midpoint + 0.5 * ds_new * RMT[DIM][i]
    end
    return midpoint_new
end
"""
$(SIGNATURES)
"""
df_refine(DIM::Integer,::AbstractVector, ::AbstractMatrix, df::AbstractVector) = ones(2^DIM) * df'
"""
$(TYPEDSIGNATURES)
"""
function vs_coarsen!(va_data::Velocity_Adaptive_Data,ka::KA{DIM,NDF})where{DIM,NDF}
    trees = ka.kdata.field.trees;kinfo = ka.kinfo
    ds = [(kinfo.config.quadrature[2*i] - kinfo.config.quadrature[2*i-1]) /
    kinfo.config.vs_trees_num[i] for i in 1:DIM]
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
            isa(ps_data,InsideSolidData) && continue
            vs_data = ps_data.vs_data
            U = @view(ps_data.prim[2:1+DIM])
            n = vs_data.vs_num
            s1 = 0.0; s2 = 0.0
            if lohner
                build_vs_index!(vsidx, vs_data, vmin, ds0, vstn, maxlevel)
                s1, s2 = vs_lohner_scales(vs_data)
            end
            resize!(coarsen_ok, n)
            @inbounds for c in 1:n
                _criterion_cell!(cdf_i, mid_i, vs_data, ps_data.ds, c)
                coarsen_ok[c] = lohner ?
                    (vs_lohner_indicator(vs_data, vsidx, c, s1, s2) < τc &&
                     local_contribution_coarsen_flag(ps_data.w, U, mid_i, cdf_i, vs_data.weight[c], kinfo)) :
                    (local_contribution_coarsen_flag(ps_data.w, U, mid_i, cdf_i, vs_data.weight[c], kinfo) &&
                     global_contribution_coarsen_flag(U, mid_i, cdf_i, vs_data.weight[c], va_data.vr, kinfo))
            end
            coarsen_grid_stream!(vs_data, coarsen_ok, ds, maxlevel) && (va_flags[id] = true)
        end
    end
    return nothing
end

function midpoint_coarsen(DIM::Integer,midpoint::AbstractVector, level::Int8, ds::AbstractVector)
    ds_new = ds / 2^level
    @. midpoint - 0.5 * ds_new * RMT[DIM][1]
end
"""
$(SIGNATURES)
"""
function df_coarsen(DIM::Integer,NDF::Integer,df::AbstractMatrix)
    df_new = zeros(NDF)
    for i = 1:2^DIM
        df_new += @view(df[i, :])
    end
    df_new /= 2^DIM
    return df_new
end


"""
$(SIGNATURES)
"""
function level_coarsen_replace!(DIM::Integer,level::AbstractVector, index::Int)
    for _ = 1:2^DIM-1
        deleteat!(level, index)
    end
    level[index] -= 1
end
"""
$(SIGNATURES)
"""
function weight_coarsen_replace!(DIM::Integer,weight::AbstractVector, index::Int)
    for _ = 1:2^DIM-1
        deleteat!(weight, index)
    end
    weight[index] *= 2^DIM
end
"""
$(SIGNATURES)
"""
function midpoint_coarsen_replace!(
    DIM::Integer,
    lnmidpoint::AbstractVector,
    midpoint_new::AbstractVector,
    vs_num::Int,
    index::Int,
)
for i = 1:DIM
        for _ = 1:2^DIM-1
            deleteat!(lnmidpoint, (i - 1) * (vs_num) + index)
        end
        lnmidpoint[(i-1)*(vs_num)+index] = midpoint_new[i]
    end
end
"""
$(SIGNATURES)
"""
function df_coarsen_replace!(DIM::Integer,NDF::Integer,lndf::AbstractVector, df_new::AbstractVector, vs_num::Int, index::Int)
    for i = 1:NDF
        for _ = 1:2^DIM-1
            deleteat!(lndf, (i - 1) * (vs_num) + index)
        end
        lndf[(i-1)*(vs_num)+index] = df_new[i]
    end
end
"""
$(SIGNATURES)
"""
function sdf_coarsen_replace!(DIM::Integer,NDF::Integer,lnsdf::AbstractVector)
    deleteat!(lnsdf, 1:(2^DIM-1)*NDF*DIM)
end
"""
$(SIGNATURES)
"""
function flux_coarsen_replace!(DIM::Integer,NDF::Integer,lnflux::AbstractVector)
    deleteat!(lnflux, 1:(2^DIM-1)*NDF)
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
            conserved_I_porjection!(vs_data,ps_data.w)
        end
    end
end

function vs_resolution(ka)
    trees = ka.kdata.field.trees
    vs_resolution(trees,ka.kinfo)
end
function vs_resolution(trees::PsTrees,kinfo)
    density_res = 0.;energy_res = 0.
    @inbounds for tree in trees.data
        for ps_data in tree
            (isa(ps_data,InsideSolidData)||ps_data.bound_enc<0)&&continue
            density_res_i,energy_res_i = vs_resolution(ps_data,kinfo)
            density_res = max(density_res_i,density_res)
            energy_res = max(energy_res_i,energy_res)
        end
    end
    density_res = MPI.Allreduce(density_res, MPI.MAX, MPI.COMM_WORLD)
    energy_res = MPI.Allreduce(energy_res, MPI.MAX, MPI.COMM_WORLD)
    return Velocity_Resolution(density_res,energy_res)
end
function vs_resolution(ps_data::PsData{DIM,NDF},kinfo) where{DIM,NDF}
    vs_data = ps_data.vs_data
    U = ps_data.prim[2:DIM+1]
    density_max = maximum(vs_data.df)
    c2 = @views [sum((U-u).^2) for u in eachrow(vs_data.midpoint)]
    energy_max = NDF==1 ? 0.5*maximum(vs_data.df.*c2) : 0.5*maximum(@view(vs_data.df[:,1]).*c2+@view(vs_data.df[:,2]))
    vs_trees_num = kinfo.config.vs_trees_num
    quadrature = kinfo.config.quadrature
    du = [quadrature[2*i]-quadrature[2*i-1] for i in 1:DIM]
    weight = reduce(*,du)/reduce(*,vs_trees_num)/2^(DIM*(kinfo.config.solver.AMR_VS_MAXLEVEL))
    return density_max*weight,energy_max*weight
end
"""
$(TYPEDSIGNATURES)
"""
function vs_adaptive_mesh_refinement!(ka;vs_balance = false)
    vr = vs_resolution(ka)
    fp = PointerWrapper(ka.kinfo.forest.p4est)
    va_flags = zeros(Bool,fp.local_num_quadrants[])
    va_data = Velocity_Adaptive_Data(vr,va_flags)
    vs_refine!(va_data,ka)
    vs_coarsen!(va_data,ka)
    changed = any(va_flags)
    if vs_balance
        changed |= vs_balance!(ka)
    end
    vs_conserved_correction!(va_data,ka)
    return Bool(MPI.Allreduce(Int(changed), +, MPI.COMM_WORLD) > 0)
end

function initial_vs_adaptive_mesh_refinement!(prim,vs_data,kinfo::KInfo{DIM,NDF}) where{DIM,NDF}
    ds = [(kinfo.config.quadrature[2*i] - kinfo.config.quadrature[2*i-1]) /
        kinfo.config.vs_trees_num[i] for i in 1:DIM]
    maxlevel = kinfo.config.solver.AMR_VS_MAXLEVEL
    ddus = [ds ./ 2.0^i for i in 0:maxlevel]      # cell size per level
    U = @view(prim[2:1+DIM])
    n = vs_data.vs_num
    refine_flags = Vector{Bool}(undef, n)
    I0buf = Vector{Float64}(undef, DIM)
    I2buf = Vector{Float64}(undef, DIM)
    @inbounds for c in 1:n
        L = vs_data.level[c]
        mid = @view(vs_data.midpoint[c, :])
        refine_flags[c] = L < maxlevel &&
            maxwellian_refine_flag(mid, ddus[L+1], U, prim, kinfo, I0buf, I2buf)
    end
    refine_grid_stream!(vs_data, refine_flags, ds)
    return nothing
end

