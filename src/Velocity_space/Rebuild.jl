# Streaming rebuild of a velocity grid under refine/coarsen.
#
# The original `vs_refine!`/`vs_coarsen!` mutate the linearized grid arrays with positional
# `deleteat!`/`insert!`, which is O(N) per refined/coarsened cell and therefore O(N^2) over a
# pass.  The streaming variants here run the decision once into a per-cell flag array and then
# emit the new grid in a single O(N) pass into fresh result arrays (kept cells copied, refined
# cells expanded, coarsened groups merged).  `sdf`/`flux` are reset to zero of the new size —
# exactly as the original (which appends zeros) — since both are recomputed by `slope!`/`flux!`
# before use.
#
# Decisions are pure functions of the *original* (pre-mutation) cell state, so precomputing
# them up front is equivalent to the original's inline evaluation.  Both an in-place reference
# (`*_grid_inplace!`, faithful extraction of the original surgery) and the streaming version
# (`*_grid_stream!`) consume the same flags; `test_vs_rebuild.jl` asserts they produce
# identical grids.

@inline function _copy_cell!(dlevel, dweight, dmid, ddf, w::Int, vs_data::AbstractVsData{DIM,NDF}, src::Int) where {DIM,NDF}
    @inbounds begin
        dlevel[w] = vs_data.level[src]
        dweight[w] = vs_data.weight[src]
        for d in 1:DIM
            dmid[w, d] = vs_data.midpoint[src, d]
        end
        for k in 1:NDF
            ddf[w, k] = vs_data.df[src, k]
        end
    end
end

@inline function _all_same_level(vs_data::AbstractVsData, index::Int, fl, nc::Int)
    @inbounds for g in 1:nc-1
        vs_data.level[index+g] != fl && return false
    end
    return true
end

# ---------------------------------------------------------------------------- refine ------

"""
$(TYPEDSIGNATURES)
In-place reference refinement (faithful extraction of the original `deleteat!`/`insert!`
surgery).  `refine_flags[i]` is the final decision for original cell `i` (already including
the `level < maxlevel` guard).  Used as the equivalence baseline for [`refine_grid_stream!`](@ref).
"""
function refine_grid_inplace!(vs_data::VsData{DIM,NDF}, refine_flags::AbstractVector{Bool}, ds) where {DIM,NDF}
    lnmidpoint = reshape(vs_data.midpoint, :)
    lndf = reshape(vs_data.df, :)
    lnsdf = reshape(vs_data.sdf, :)
    lnflux = reshape(vs_data.flux, :)
    midpoint_index = Vector{Int}(undef, DIM)
    df_index = Vector{Int}(undef, NDF)
    changed = false
    index = 1; oi = 0
    while index < vs_data.vs_num + 1
        oi += 1
        @inbounds for i in 1:DIM
            midpoint_index[i] = (i - 1) * vs_data.vs_num + index
        end
        @inbounds for i in 1:NDF
            df_index[i] = (i - 1) * vs_data.vs_num + index
        end
        midpoint = @view(lnmidpoint[midpoint_index])
        df = @view(lndf[df_index])
        if refine_flags[oi]
            midpoint_new = midpoint_refine(DIM, midpoint, vs_data.level[index], ds)
            df_new = df_refine(DIM, midpoint, midpoint_new, df)
            vs_data.vs_num += 2^DIM - 1
            level_refine_replace!(DIM, vs_data.level, index)
            weight_refine_replace!(DIM, vs_data.weight, index)
            midpoint_refine_replace!(DIM, lnmidpoint, midpoint_new, vs_data.vs_num, index)
            df_refine_replace!(DIM, NDF, lndf, df_new, vs_data.vs_num, index)
            sdf_refine_replace!(DIM, NDF, lnsdf)
            flux_refine_replace!(DIM, NDF, lnflux)
            index += 2^DIM - 1
            changed = true
        end
        index += 1
    end
    vs_data.midpoint = reshape(lnmidpoint, vs_data.vs_num, DIM)
    vs_data.df = reshape(lndf, vs_data.vs_num, NDF)
    vs_data.sdf = reshape(lnsdf, vs_data.vs_num, NDF, DIM)
    vs_data.flux = reshape(lnflux, vs_data.vs_num, NDF)
    return changed
end

"""
$(TYPEDSIGNATURES)
Streaming refinement: single O(N) pass into fresh arrays.  `refine_flags[i]` decides original
cell `i`.  Children are emitted in `RMT` order with the parent's `df` injected (matching
`df_refine`), parent weight split by `2^DIM`, and `sdf`/`flux` reset to zero.
"""
function refine_grid_stream!(vs_data::VsData{DIM,NDF}, refine_flags::AbstractVector{Bool}, ds) where {DIM,NDF}
    n = vs_data.vs_num
    nc = 2^DIM
    nnew = n
    @inbounds for i in 1:n
        refine_flags[i] && (nnew += nc - 1)
    end
    nnew == n && return false
    dlevel = Vector{Int8}(undef, nnew)
    dweight = Vector{Float64}(undef, nnew)
    dmid = Matrix{Float64}(undef, nnew, DIM)
    ddf = Matrix{Float64}(undef, nnew, NDF)
    w = 0
    @inbounds for i in 1:n
        if refine_flags[i]
            L = vs_data.level[i]
            wi = vs_data.weight[i] / nc
            for c in 1:nc
                w += 1
                dlevel[w] = L + 1
                dweight[w] = wi
                for d in 1:DIM
                    dmid[w, d] = vs_data.midpoint[i, d] + 0.5 * (ds[d] / 2^(L + 1)) * RMT[DIM][c][d]
                end
                for k in 1:NDF
                    ddf[w, k] = vs_data.df[i, k]
                end
            end
        else
            w += 1
            _copy_cell!(dlevel, dweight, dmid, ddf, w, vs_data, i)
        end
    end
    vs_data.level = dlevel
    vs_data.weight = dweight
    vs_data.midpoint = dmid
    vs_data.df = ddf
    vs_data.sdf = zeros(nnew, NDF, DIM)
    vs_data.flux = zeros(nnew, NDF)
    vs_data.vs_num = nnew
    return true
end

# --------------------------------------------------------------------------- coarsen ------

"""
$(TYPEDSIGNATURES)
In-place reference coarsening (faithful extraction of the original surgery).  `coarsen_ok[i]`
is the per-original-cell willingness to coarsen; a `2^DIM`-aligned same-level sibling group is
merged iff every member is willing.  `maxlevel` sizes the alignment bookkeeping `flag`.
"""
function coarsen_grid_inplace!(vs_data::VsData{DIM,NDF}, coarsen_ok::AbstractVector{Bool}, ds, maxlevel::Integer) where {DIM,NDF}
    lnmidpoint = reshape(vs_data.midpoint, :)
    lndf = reshape(vs_data.df, :)
    lnsdf = reshape(vs_data.sdf, :)
    lnflux = reshape(vs_data.flux, :)
    midpoint_index = Matrix{Int}(undef, 2^DIM, DIM)
    df_index = Matrix{Int}(undef, 2^DIM, NDF)
    flag = zeros(maxlevel)
    nc = 2^DIM
    changed = false
    index = 1; op = 1
    while index < vs_data.vs_num + 1
        first_level = vs_data.level[index]
        if first_level > 0
            if flag[first_level] % 1 == 0.0 && index + nc - 1 <= vs_data.vs_num &&
               _all_same_level(vs_data, index, first_level, nc)
                grp = true
                @inbounds for g in 0:nc-1
                    coarsen_ok[op+g] || (grp = false; break)
                end
                if grp
                    for i in axes(midpoint_index, 2)
                        midpoint_index[:, i] .=
                            (i-1)*vs_data.vs_num+index:(i-1)*vs_data.vs_num+index+nc-1
                    end
                    for i in axes(df_index, 2)
                        df_index[:, i] .=
                            (i-1)*vs_data.vs_num+index:(i-1)*vs_data.vs_num+index+nc-1
                    end
                    midpoint = @view(lnmidpoint[midpoint_index])
                    df = @view(lndf[df_index])
                    midpoint_new = midpoint_coarsen(DIM, @view(midpoint[1, :]), first_level, ds)
                    df_new = df_coarsen(DIM, NDF, df)
                    vs_data.vs_num -= nc - 1
                    level_coarsen_replace!(DIM, vs_data.level, index)
                    weight_coarsen_replace!(DIM, vs_data.weight, index)
                    midpoint_coarsen_replace!(DIM, lnmidpoint, midpoint_new, vs_data.vs_num, index)
                    df_coarsen_replace!(DIM, NDF, lndf, df_new, vs_data.vs_num, index)
                    sdf_coarsen_replace!(DIM, NDF, lnsdf)
                    flux_coarsen_replace!(DIM, NDF, lnflux)
                    changed = true
                else
                    index += nc - 1
                end
                op += nc
                if first_level > 1
                    for i in 1:first_level-1
                        flag[i] += 1 / 2^(DIM * (first_level - i))
                    end
                end
            else
                for i in 1:first_level
                    flag[i] += 1 / 2^(DIM * (first_level - i + 1))
                end
                op += 1
            end
        else
            op += 1
        end
        index += 1
    end
    vs_data.sdf = reshape(lnsdf, vs_data.vs_num, NDF, DIM)
    vs_data.flux = reshape(lnflux, vs_data.vs_num, NDF)
    vs_data.midpoint = reshape(lnmidpoint, vs_data.vs_num, DIM)
    vs_data.df = reshape(lndf, vs_data.vs_num, NDF)
    return changed
end

"""
$(TYPEDSIGNATURES)
Streaming coarsening: single read-only pass over the source grid driving the same alignment
walk as the original, emitting merged/kept cells into fresh arrays.  Equivalent to
[`coarsen_grid_inplace!`](@ref).
"""
function coarsen_grid_stream!(vs_data::VsData{DIM,NDF}, coarsen_ok::AbstractVector{Bool}, ds, maxlevel::Integer) where {DIM,NDF}
    n = vs_data.vs_num
    nc = 2^DIM
    dlevel = Vector{Int8}(undef, n)
    dweight = Vector{Float64}(undef, n)
    dmid = Matrix{Float64}(undef, n, DIM)
    ddf = Matrix{Float64}(undef, n, NDF)
    flag = zeros(maxlevel)
    changed = false
    w = 0; index = 1
    @inbounds while index <= n
        fl = vs_data.level[index]
        if fl > 0
            if flag[fl] % 1 == 0.0 && index + nc - 1 <= n && _all_same_level(vs_data, index, fl, nc)
                grp = true
                for g in 0:nc-1
                    coarsen_ok[index+g] || (grp = false; break)
                end
                if grp
                    w += 1
                    dlevel[w] = fl - 1
                    dweight[w] = vs_data.weight[index] * nc
                    for d in 1:DIM
                        dmid[w, d] = vs_data.midpoint[index, d] - 0.5 * (ds[d] / 2^fl) * RMT[DIM][1][d]
                    end
                    for k in 1:NDF
                        s = 0.0
                        for g in 0:nc-1
                            s += vs_data.df[index+g, k]
                        end
                        ddf[w, k] = s / nc
                    end
                    changed = true
                else
                    for g in 0:nc-1
                        w += 1
                        _copy_cell!(dlevel, dweight, dmid, ddf, w, vs_data, index + g)
                    end
                end
                index += nc
                if fl > 1
                    for i in 1:fl-1
                        flag[i] += 1 / 2^(DIM * (fl - i))
                    end
                end
            else
                w += 1
                _copy_cell!(dlevel, dweight, dmid, ddf, w, vs_data, index)
                for i in 1:fl
                    flag[i] += 1 / 2^(DIM * (fl - i + 1))
                end
                index += 1
            end
        else
            w += 1
            _copy_cell!(dlevel, dweight, dmid, ddf, w, vs_data, index)
            index += 1
        end
    end
    vs_data.level = resize!(dlevel, w)
    vs_data.weight = resize!(dweight, w)
    vs_data.midpoint = dmid[1:w, :]
    vs_data.df = ddf[1:w, :]
    vs_data.sdf = zeros(w, NDF, DIM)
    vs_data.flux = zeros(w, NDF)
    vs_data.vs_num = w
    return changed
end
