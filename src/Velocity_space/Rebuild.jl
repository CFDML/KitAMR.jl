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

@inline function _push_unique_vs_neighbor!(buf::Vector{Int}, nb::Integer, center::Integer)
    nb == 0 && return nothing
    nb == center && return nothing
    @inbounds for x in buf
        x == nb && return nothing
    end
    push!(buf, Int(nb))
    return nothing
end

function _linear_refine_gradient!(
    grad::AbstractMatrix{Float64},
    vs_data::AbstractVsData{DIM,NDF},
    idx::VsNeighborIndex{DIM},
    src::Integer,
    ds,
    neighbors::Vector{Int},
    normal::AbstractMatrix{Float64},
    rhs::AbstractVector{Float64},
    xbuf::AbstractVector{Float64},
) where {DIM,NDF}
    fill!(grad, 0.0)
    empty!(neighbors)
    @inbounds for d in 1:DIM
        _push_unique_vs_neighbor!(neighbors, vs_face_neighbor(idx, vs_data, src, d, -1), src)
        _push_unique_vs_neighbor!(neighbors, vs_face_neighbor(idx, vs_data, src, d, 1), src)
    end
    isempty(neighbors) && return grad

    fill!(normal, 0.0)
    @inbounds for nb in neighbors
        dist2 = 0.0
        for d in 1:DIM
            h = ds[d] / 2.0^vs_data.level[src]
            xbuf[d] = (vs_data.midpoint[nb, d] - vs_data.midpoint[src, d]) / h
            dist2 += xbuf[d]^2
        end
        dist2 <= EPS && continue
        wt = 1.0 / (1.0 + dist2)
        for a in 1:DIM
            xa = xbuf[a]
            for b in a:DIM
                normal[a, b] += wt * xa * xbuf[b]
            end
        end
    end
    @inbounds for a in 1:DIM
        for b in 1:a-1
            normal[a, b] = normal[b, a]
        end
    end
    traceN = 0.0
    @inbounds for d in 1:DIM
        traceN += normal[d, d]
    end
    traceN <= EPS && return grad
    ridge = 1e-12 * max(traceN, 1.0)
    @inbounds for d in 1:DIM
        normal[d, d] += ridge
    end

    @inbounds for k in 1:NDF
        fill!(rhs, 0.0)
        f0 = vs_data.df[src, k]
        for nb in neighbors
            dist2 = 0.0
            for d in 1:DIM
                h = ds[d] / 2.0^vs_data.level[src]
                xbuf[d] = (vs_data.midpoint[nb, d] - vs_data.midpoint[src, d]) / h
                dist2 += xbuf[d]^2
            end
            dist2 <= EPS && continue
            wt = 1.0 / (1.0 + dist2)
            df = vs_data.df[nb, k] - f0
            for a in 1:DIM
                rhs[a] += wt * xbuf[a] * df
            end
        end
        slope = try
            Symmetric(normal) \ rhs
        catch
            nothing
        end
        slope === nothing && continue
        for d in 1:DIM
            grad[k, d] = slope[d]
        end
    end
    return grad
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

"""
$(TYPEDSIGNATURES)
Streaming refinement with neighbor-based linear prolongation.  The refinement decisions and
the neighbor index are both defined on the original grid.  For each refined parent, a small
least-squares affine model is fitted from its face neighbors in scaled local coordinates,
then evaluated at the child centers.  No limiter is applied here; conservative
I-projection/shaving can be applied afterwards by the caller.
"""
function refine_grid_stream_linear!(
    vs_data::VsData{DIM,NDF},
    refine_flags::AbstractVector{Bool},
    ds,
    idx::VsNeighborIndex{DIM},
) where {DIM,NDF}
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
    neighbors = Int[]
    normal = zeros(Float64, DIM, DIM)
    rhs = zeros(Float64, DIM)
    xbuf = zeros(Float64, DIM)
    grad = zeros(Float64, NDF, DIM)
    w = 0
    @inbounds for i in 1:n
        if refine_flags[i]
            L = vs_data.level[i]
            wi = vs_data.weight[i] / nc
            _linear_refine_gradient!(grad, vs_data, idx, i, ds, neighbors, normal, rhs, xbuf)
            for c in 1:nc
                w += 1
                dlevel[w] = L + 1
                dweight[w] = wi
                for d in 1:DIM
                    dmid[w, d] = vs_data.midpoint[i, d] + 0.5 * (ds[d] / 2^(L + 1)) * RMT[DIM][c][d]
                end
                for k in 1:NDF
                    val = vs_data.df[i, k]
                    for d in 1:DIM
                        h = ds[d] / 2.0^L
                        val += grad[k, d] * ((dmid[w, d] - vs_data.midpoint[i, d]) / h)
                    end
                    ddf[w, k] = val
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
