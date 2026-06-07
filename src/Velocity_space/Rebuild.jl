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
# them up front is equivalent to the original's inline evaluation.  Equivalence with the
# original positional-surgery refine/coarsen was verified bit-for-bit before that code was
# retired.

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
Streaming refinement: single O(N) pass into fresh arrays.  `refine_flags[i]` decides original
cell `i`.  Children are emitted in `RMT` order with the parent's `df` injected, parent weight
split by `2^DIM`, and `sdf`/`flux` reset to zero.
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
Streaming coarsening: single read-only pass over the source grid driving the alignment walk,
emitting merged/kept cells into fresh arrays.  `coarsen_ok[i]` is the per-cell willingness to
coarsen; a `2^DIM`-aligned same-level sibling group is merged iff every member is willing.
`maxlevel` sizes the alignment bookkeeping `flag`.
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
