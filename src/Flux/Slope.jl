"""
$(TYPEDSIGNATURES)
"""
function vanleer(sL::AbstractArray, sR::AbstractArray)
    SL = abs.(sL)
    SR = abs.(sR)
    @. (sign(sL) + sign(sR)) * SL * SR / (SL + SR + EPS)
end
"""
$(TYPEDSIGNATURES)
"""
function vanleer(sL::Real, sR::Real)
    SL = abs(sL)
    SR = abs(sR)
    (sign(sL) + sign(sR)) * SL * SR / (SL + SR + EPS)
end
"""
$(TYPEDSIGNATURES)
"""
function minmod(sL::Real,sR::Real)
    SL = abs(sL)
    SR = abs(sR)
    0.5(sign(sL)+sign(sR))*min(SL,SR)
end
"""
$(TYPEDSIGNATURES)
Difference through all velocity cells in two neighboring physical cells.
"""
function diff_vs!(vs_data::AbstractVsData{DIM,NDF}, vs_data_n::AbstractVsData{DIM,NDF}, dsL::Float64, sL::AbstractMatrix) where{DIM,NDF}
    index = 1
    flag = 0.0
    level = vs_data.level
    level_n = vs_data_n.level
    df = vs_data.df
    dfn = vs_data_n.df
    @inbounds for i = 1:vs_data.vs_num
        if level[i] == level_n[index]
            for j = 1:NDF
                sL[i, j] += (df[i, j] - dfn[index, j]) / dsL
            end
            index += 1
        elseif level[i] < level_n[index]
            while flag != 1.0
                for j = 1:NDF
                    sL[i, j] +=
                        (df[i, j] - dfn[index, j]) / 2^(DIM * (level_n[index] - level[i])) /
                        dsL
                end
                flag += 1 / 2^(DIM * (level_n[index] - level[i]))
                index += 1
            end
            flag = 0.0
        else
            for j = 1:NDF
                sL[i, j] += (df[i, j] - dfn[index, j]) / dsL
            end
            flag += 1 / 2^(DIM * (level[i] - level_n[index]))
            if flag == 1.0
                index += 1
                flag = 0.0
            end
        end
    end
end
"""
$(TYPEDSIGNATURES)
"""
function update_slope_bound_vs!(
    vs_data::AbstractVsData{DIM,NDF},
    neighbor_datas::AbstractVector,
    ds::Float64,
    dir::Int,
    ws_sdf::Matrix{Float64},
) where{DIM,NDF}
    vsn = vs_data.vs_num
    sdf = @view ws_sdf[1:vsn, :]
    fill!(sdf, 0.0)
    N = length(neighbor_datas)
    for j = 1:N
        neighbor_data = neighbor_datas[j].vs_data
        diff_vs!(vs_data, neighbor_data, ds, sdf)
    end
    sdf ./= N
    vs_data.sdf[:, :, dir] .= sdf
    return nothing
end
"""
$(TYPEDSIGNATURES)
"""
function update_slope_inner_vs!(
    vs_data::AbstractVsData{DIM,NDF},
    Ldata::Vector,
    Rdata::Vector,
    dsL::Float64,
    dsR::Float64,
    dir::Int,
    ws_sL::Matrix{Float64},
    ws_sR::Matrix{Float64},
) where{DIM,NDF}
    vsn = vs_data.vs_num
    sL = @view ws_sL[1:vsn, :]
    sR = @view ws_sR[1:vsn, :]
    fill!(sL, 0.0); fill!(sR, 0.0)
    nL = length(Ldata); nR = length(Rdata)
    for j in 1:nL
        L_data = Ldata[j].vs_data
        diff_vs!(vs_data, L_data, dsL, sL)
    end
    for j in 1:nR
        R_data = Rdata[j].vs_data
        diff_vs!(vs_data, R_data, dsR, sR)
    end
    sL ./= nL; sR ./= nR
    @. vs_data.sdf[:, :, dir] = minmod(sL, sR)
    return nothing
end
"""
$(TYPEDSIGNATURES)
"""
function update_slope_inner_ps!(
    ps_data::PsData{DIM,NDF},
    Ldata::AbstractVector,
    Rdata::AbstractVector,
    dsL::Float64,
    dsR::Float64,
    dir::Int,
    ws_swL::Vector{Float64},
    ws_swR::Vector{Float64},
) where{DIM,NDF}
    fill!(ws_swL, 0.0); fill!(ws_swR, 0.0)
    nL = length(Ldata); nR = length(Rdata)
    for j in 1:nL
        @. ws_swL += (ps_data.w - Ldata[j].w) / dsL
    end
    for j in 1:nR
        @. ws_swR += (ps_data.w - Rdata[j].w) / dsR
    end
    ws_swL ./= nL; ws_swR ./= nR
    for j in eachindex(ws_swL)
        ps_data.sw[j, dir] = minmod(ws_swL[j], ws_swR[j])
    end
    return nothing
end

"""
$(TYPEDSIGNATURES)
"""
function update_slope_bound_ps!(
    ps_data::PsData{DIM,NDF},
    Ldata::AbstractVector,
    dsL::Float64,
    dir::Int,
    ws_swL::Vector{Float64},
) where{DIM,NDF}
    fill!(ws_swL, 0.0)
    nL = length(Ldata)
    for j in 1:nL
        @. ws_swL += (ps_data.w - Ldata[j].w) / dsL
    end
    ws_swL ./= nL
    ps_data.sw[:, dir] .= ws_swL
    return nothing
end

# ---------------------------------------------------------------------------
# Transverse-projection variants used in the per-level correction pass of
# `slope!`. For a fine cell adjacent to a coarse neighbour, the coarse
# neighbour's midpoint is offset from the fine cell's midpoint in the
# direction(s) perpendicular to `dir` (the face normal). The raw differences
# `(ps_data.w - n.w)/ds` therefore include a spurious transverse-gradient
# contribution `≈ (Δm_t) · ∂w/∂t / ds`. We remove it by projecting the
# neighbour's `w` onto the fine cell's transverse coordinate using the
# neighbour's own transverse slope `sw[:, t]` (assumed already final, i.e.
# `slope!` updates levels coarsest-to-finest).
# ---------------------------------------------------------------------------

@inline function _ps_has_transverse_offset(ps_data::PsData{DIM}, neighbours::AbstractVector, dir::Int) where{DIM}
    DIM == 1 && return false
    n_count = length(neighbours)
    n_count == 0 && return false
    # Compare the *averaged* transverse midpoint of the neighbour set with
    # `ps_data.midpoint`. For a single coarse neighbour the average equals
    # its midpoint and may differ from `ps_data` (the offset we want to
    # correct). For a hanging-fine neighbour list the fine subcells'
    # midpoints are symmetric around `ps_data.midpoint` in the transverse
    # direction(s), so the average matches `ps_data` and no projection is
    # needed (and would otherwise rely on the fine neighbours' own sw,
    # which has not yet been computed in the coarsest-to-finest sweep).
    @inbounds for t in 1:DIM
        t == dir && continue
        avg_t = 0.0
        for j in 1:n_count
            avg_t += neighbours[j].midpoint[t]
        end
        avg_t /= n_count
        avg_t != ps_data.midpoint[t] && return true
    end
    return false
end

"""
$(TYPEDSIGNATURES)
Returns `true` if any neighbour in `Ldata` or `Rdata` has a transverse-midpoint
offset relative to `ps_data` in `dir`. Used to decide whether the transverse
projection is needed at all on a given face (1-D fine cells skip the work).
"""
@inline function needs_transverse_correction(ps_data::PsData, Ldata::AbstractVector, Rdata::AbstractVector, dir::Int)
    _ps_has_transverse_offset(ps_data, Ldata, dir) ||
    _ps_has_transverse_offset(ps_data, Rdata, dir)
end

"""
$(TYPEDSIGNATURES)
Inner PS-slope reconstruction with per-side transverse projection. The two
sides are treated independently: each side projects only if its neighbour
set is geometrically offset from `ps_data` in the transverse direction
(checked via [`_ps_has_transverse_offset`](@ref)). This skips hanging-fine
sides — those average to `ps_data.midpoint` and would otherwise depend on
fine neighbours' sw that has not yet been computed in the coarsest-to-finest
sweep. The vanleer limiter on the final L/R slopes is preserved.
"""
function update_slope_inner_ps_transverse!(
    ps_data::PsData{DIM,NDF},
    Ldata::AbstractVector,
    Rdata::AbstractVector,
    dsL::Float64,
    dsR::Float64,
    dir::Int,
    ws_swL::Vector{Float64},
    ws_swR::Vector{Float64},
) where{DIM,NDF}
    fill!(ws_swL, 0.0); fill!(ws_swR, 0.0)
    nL = length(Ldata); nR = length(Rdata)
    projL = _ps_has_transverse_offset(ps_data, Ldata, dir)
    projR = _ps_has_transverse_offset(ps_data, Rdata, dir)
    @inbounds for j in 1:nL
        n = Ldata[j]
        for i in eachindex(ws_swL)
            w_proj = n.w[i]
            if projL
                for t in 1:DIM
                    t == dir && continue
                    w_proj += (ps_data.midpoint[t] - n.midpoint[t]) * n.sw[i, t]
                end
            end
            ws_swL[i] += (ps_data.w[i] - w_proj) / dsL
        end
    end
    @inbounds for j in 1:nR
        n = Rdata[j]
        for i in eachindex(ws_swR)
            w_proj = n.w[i]
            if projR
                for t in 1:DIM
                    t == dir && continue
                    w_proj += (ps_data.midpoint[t] - n.midpoint[t]) * n.sw[i, t]
                end
            end
            ws_swR[i] += (ps_data.w[i] - w_proj) / dsR
        end
    end
    ws_swL ./= nL; ws_swR ./= nR
    @inbounds for j in eachindex(ws_swL)
        ps_data.sw[j, dir] = minmod(ws_swL[j], ws_swR[j])
    end
    return nothing
end

"""
$(TYPEDSIGNATURES)
Transverse-projected difference through all velocity cells in two neighbouring
physical cells. Variant of [`diff_vs!`](@ref) that subtracts the neighbour's
transverse df-slope `sdf[:,:,t]` to remove the spurious contribution from the
transverse midpoint offset between `ps_data` and `vs_data_n`.

`dm_t` carries the transverse midpoint offsets (length DIM, with `dm_t[dir]=0`).
"""
function diff_vs_transverse!(
    vs_data::AbstractVsData{DIM,NDF},
    vs_data_n::AbstractVsData{DIM,NDF},
    dsL::Float64,
    dm_t::AbstractVector{Float64},
    sL::AbstractMatrix,
) where{DIM,NDF}
    index = 1
    flag = 0.0
    level = vs_data.level
    level_n = vs_data_n.level
    df = vs_data.df
    dfn = vs_data_n.df
    sdfn = vs_data_n.sdf
    @inbounds for i = 1:vs_data.vs_num
        if level[i] == level_n[index]
            for j = 1:NDF
                proj = dfn[index, j]
                for t in 1:DIM
                    proj += dm_t[t] * sdfn[index, j, t]
                end
                sL[i, j] += (df[i, j] - proj) / dsL
            end
            index += 1
        elseif level[i] < level_n[index]
            while flag != 1.0
                for j = 1:NDF
                    proj = dfn[index, j]
                    for t in 1:DIM
                        proj += dm_t[t] * sdfn[index, j, t]
                    end
                    sL[i, j] +=
                        (df[i, j] - proj) / 2^(DIM * (level_n[index] - level[i])) /
                        dsL
                end
                flag += 1 / 2^(DIM * (level_n[index] - level[i]))
                index += 1
            end
            flag = 0.0
        else
            for j = 1:NDF
                proj = dfn[index, j]
                for t in 1:DIM
                    proj += dm_t[t] * sdfn[index, j, t]
                end
                sL[i, j] += (df[i, j] - proj) / dsL
            end
            flag += 1 / 2^(DIM * (level[i] - level_n[index]))
            if flag == 1.0
                index += 1
                flag = 0.0
            end
        end
    end
end

"""
$(TYPEDSIGNATURES)
Inner VS-slope reconstruction with transverse projection. Counterpart of
[`update_slope_inner_ps_transverse!`](@ref) for `vs_data.sdf[:, :, dir]`.
"""
function update_slope_inner_vs_transverse!(
    ps_data::PsData{DIM,NDF},
    Ldata::Vector,
    Rdata::Vector,
    dsL::Float64,
    dsR::Float64,
    dir::Int,
    ws_sL::Matrix{Float64},
    ws_sR::Matrix{Float64},
) where{DIM,NDF}
    vs_data = ps_data.vs_data
    vsn = vs_data.vs_num
    sL = @view ws_sL[1:vsn, :]
    sR = @view ws_sR[1:vsn, :]
    fill!(sL, 0.0); fill!(sR, 0.0)
    nL = length(Ldata); nR = length(Rdata)
    projL = _ps_has_transverse_offset(ps_data, Ldata, dir)
    projR = _ps_has_transverse_offset(ps_data, Rdata, dir)
    dm = zeros(Float64, DIM)
    @inbounds for j in 1:nL
        n = Ldata[j]
        if projL
            for t in 1:DIM
                dm[t] = (t == dir) ? 0.0 : (ps_data.midpoint[t] - n.midpoint[t])
            end
            diff_vs_transverse!(vs_data, n.vs_data, dsL, dm, sL)
        else
            diff_vs!(vs_data, n.vs_data, dsL, sL)
        end
    end
    @inbounds for j in 1:nR
        n = Rdata[j]
        if projR
            for t in 1:DIM
                dm[t] = (t == dir) ? 0.0 : (ps_data.midpoint[t] - n.midpoint[t])
            end
            diff_vs_transverse!(vs_data, n.vs_data, dsR, dm, sR)
        else
            diff_vs!(vs_data, n.vs_data, dsR, sR)
        end
    end
    sL ./= nL; sR ./= nR
    @. vs_data.sdf[:, :, dir] = minmod(sL, sR)
    return nothing
end

# ---------------------------------------------------------------------------
# One-sided transverse-projection variants (used at "corner" cells where one
# face is a domain boundary and the opposite face has a coarse neighbour).
# Mirror the structure of update_slope_bound_ps!/_vs! but apply the
# transverse projection on the single available neighbour list.
# ---------------------------------------------------------------------------

"""
$(TYPEDSIGNATURES)
One-sided PS-slope reconstruction with transverse projection. Used when only
one side of `dir` has neighbours (the other side is a domain boundary) and
that side is a coarse-level neighbour. Writes `ps_data.sw[:, dir]` directly.
"""
function update_slope_bound_ps_transverse!(
    ps_data::PsData{DIM,NDF},
    Ldata::AbstractVector,
    dsL::Float64,
    dir::Int,
    ws_swL::Vector{Float64},
) where{DIM,NDF}
    fill!(ws_swL, 0.0)
    nL = length(Ldata)
    @inbounds for j in 1:nL
        n = Ldata[j]
        for i in eachindex(ws_swL)
            w_proj = n.w[i]
            for t in 1:DIM
                t == dir && continue
                w_proj += (ps_data.midpoint[t] - n.midpoint[t]) * n.sw[i, t]
            end
            ws_swL[i] += (ps_data.w[i] - w_proj) / dsL
        end
    end
    ws_swL ./= nL
    ps_data.sw[:, dir] .= ws_swL
    return nothing
end

"""
$(TYPEDSIGNATURES)
One-sided VS-slope reconstruction with transverse projection. Counterpart of
[`update_slope_bound_ps_transverse!`](@ref) for `vs_data.sdf[:, :, dir]`.
"""
function update_slope_bound_vs_transverse!(
    ps_data::PsData{DIM,NDF},
    Ldata::AbstractVector,
    dsL::Float64,
    dir::Int,
    ws_sL::Matrix{Float64},
) where{DIM,NDF}
    vs_data = ps_data.vs_data
    vsn = vs_data.vs_num
    sL = @view ws_sL[1:vsn, :]
    fill!(sL, 0.0)
    nL = length(Ldata)
    dm = zeros(Float64, DIM)
    @inbounds for j in 1:nL
        n = Ldata[j]
        for t in 1:DIM
            dm[t] = (t == dir) ? 0.0 : (ps_data.midpoint[t] - n.midpoint[t])
        end
        diff_vs_transverse!(vs_data, n.vs_data, dsL, dm, sL)
    end
    sL ./= nL
    vs_data.sdf[:, :, dir] .= sL
    return nothing
end


# Inner
"""
$(TYPEDSIGNATURES)
"""
function update_slope!(
    ::Val{1},
    ::Val{1},
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
    if Ldata[1].bound_enc<0&&Rdata[1].bound_enc<0
        ps_data.vs_data.sdf[:,:,dir] .= 0.
    elseif Ldata[1].bound_enc<0
        ds = ps_data.midpoint[dir]-Rdata[1].midpoint[dir]
        update_slope_bound_vs!(ps_data.vs_data, Rdata, ds, dir, ws_sL)
        # update_slope_bound_ps!(ps_data, Rdata, ds, dir, ws_swL)
    elseif Rdata[1].bound_enc<0
        ds = ps_data.midpoint[dir]-Ldata[1].midpoint[dir]
        update_slope_bound_vs!(ps_data.vs_data, Ldata, ds, dir, ws_sL)
        # update_slope_bound_ps!(ps_data, Ldata, ds, dir, ws_swL)
    else
        ds = ps_data.ds[dir]
        update_slope_inner_vs!(ps_data.vs_data, Ldata, Rdata, ds, -ds, dir, ws_sL, ws_sR)
        # update_slope_inner_ps!(ps_data, Ldata, Rdata, ds, -ds, dir, ws_swL, ws_swR)
    end
    return nothing
end
"""
$(TYPEDSIGNATURES)
"""
function update_slope!(
    ::Val{1},
    ::NeighborNum,
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
    update_slope_inner_vs!(ps_data.vs_data, Ldata, Rdata, ds, -0.75 * ds, dir, ws_sL, ws_sR)
    # update_slope_inner_ps!(ps_data, Ldata, Rdata, ds, -0.75 * ds, dir, ws_swL, ws_swR)
end
"""
$(TYPEDSIGNATURES)
"""
function update_slope!(
    ::NeighborNum,
    ::Val{1},
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
    update_slope_inner_vs!(ps_data.vs_data, Ldata, Rdata, 0.75 * ds, -ds, dir, ws_sL, ws_sR)
    # update_slope_inner_ps!(ps_data, Ldata, Rdata, 0.75 * ds, -ds, dir, ws_swL, ws_swR)
end
"""
$(TYPEDSIGNATURES)
"""
function update_slope!(
    ::NeighborNum,
    ::NeighborNum,
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
    update_slope_inner_vs!(ps_data.vs_data, Ldata, Rdata, 0.75 * ds, -0.75 * ds, dir, ws_sL, ws_sR)
    # update_slope_inner_ps!(ps_data, Ldata, Rdata, 0.75 * ds, -0.75 * ds, dir, ws_swL, ws_swR)
end
"""
$(TYPEDSIGNATURES)
"""
function update_slope!(
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
    update_slope_inner_vs!(ps_data.vs_data, Ldata, Rdata, 1.5 * ds, -1.5 * ds, dir, ws_sL, ws_sR)
    # update_slope_inner_ps!(ps_data, Ldata, Rdata, 1.5 * ds, -1.5 * ds, dir, ws_swL, ws_swR)
end
"""
$(TYPEDSIGNATURES)
"""
function update_slope!(
    ::Val{1},
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
    update_slope_inner_vs!(ps_data.vs_data, Ldata, Rdata, ds, -1.5 * ds, dir, ws_sL, ws_sR)
    # update_slope_inner_ps!(ps_data, Ldata, Rdata, ds, -1.5 * ds, dir, ws_swL, ws_swR)
end
"""
$(TYPEDSIGNATURES)
"""
function update_slope!(
    ::Val{-1},
    ::Val{1},
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
    update_slope_inner_vs!(ps_data.vs_data, Ldata, Rdata, 1.5 * ds, -ds, dir, ws_sL, ws_sR)
    # update_slope_inner_ps!(ps_data, Ldata, Rdata, 1.5 * ds, -ds, dir, ws_swL, ws_swR)
end
"""
$(TYPEDSIGNATURES)
"""
function update_slope!(
    ::NeighborNum,
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
    update_slope_inner_vs!(ps_data.vs_data, Ldata, Rdata, 0.75 * ds, -1.5 * ds, dir, ws_sL, ws_sR)
    # update_slope_inner_ps!(ps_data, Ldata, Rdata, 0.75 * ds, -1.5 * ds, dir, ws_swL, ws_swR)
end
"""
$(TYPEDSIGNATURES)
"""
function update_slope!(
    ::Val{-1},
    ::NeighborNum,
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
    update_slope_inner_vs!(ps_data.vs_data, Ldata, Rdata, 1.5 * ds, -0.75 * ds, dir, ws_sL, ws_sR)
    # update_slope_inner_ps!(ps_data, Ldata, Rdata, 1.5 * ds, -0.75 * ds, dir, ws_swL, ws_swR)
end

# Boundary
"""
$(TYPEDSIGNATURES)
"""
function update_slope!(
    ::Val{0},
    ::Val{1},
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
    ds = ps_data.midpoint[dir]-Rdata[1].midpoint[dir]
    vs_data = ps_data.vs_data
    update_slope_bound_vs!(vs_data, Rdata, ds, dir, ws_sL)
    # update_slope_bound_ps!(ps_data, Rdata, ds, dir, ws_swL)
end
"""
$(TYPEDSIGNATURES)
"""
function update_slope!(
    ::Val{1},
    ::Val{0},
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
    ds = ps_data.midpoint[dir]-Ldata[1].midpoint[dir]
    vs_data = ps_data.vs_data
    update_slope_bound_vs!(vs_data, Ldata, ds, dir, ws_sL)
    # update_slope_bound_ps!(ps_data, Ldata, ds, dir, ws_swL)
end
"""
$(TYPEDSIGNATURES)
"""
function update_slope!(
    ::NeighborNum,
    ::Val{0},
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
    ds = ps_data.midpoint[dir]-Ldata[1].midpoint[dir]
    update_slope_bound_vs!(ps_data.vs_data, Ldata, ds, dir, ws_sL)
    # update_slope_bound_ps!(ps_data, Ldata, ds, dir, ws_swL)
end
"""
$(TYPEDSIGNATURES)
"""
function update_slope!(
    ::Val{0},
    ::NeighborNum,
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
    ds = ps_data.midpoint[dir]-Rdata[1].midpoint[dir]
    update_slope_bound_vs!(ps_data.vs_data, Rdata, ds, dir, ws_sL)
    # update_slope_bound_ps!(ps_data, Rdata, ds, dir, ws_swL)
end
"""
$(TYPEDSIGNATURES)
"""
function update_slope!(
    ::Val{0},
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
    ds = ps_data.midpoint[dir]-Rdata[1].midpoint[dir]
    update_slope_bound_vs!(ps_data.vs_data, Rdata, ds, dir, ws_sL)
    # update_slope_bound_ps!(ps_data, Rdata, ds, dir, ws_swL)
end
"""
$(TYPEDSIGNATURES)
"""
function update_slope!(
    ::Val{-1},
    ::Val{0},
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
    ds = ps_data.midpoint[dir]-Ldata[1].midpoint[dir]
    update_slope_bound_vs!(ps_data.vs_data, Ldata, ds, dir, ws_sL)
    # update_slope_bound_ps!(ps_data, Ldata, ds, dir, ws_swL)
end

"""
$(TYPEDSIGNATURES)
"""
function update_slope!(ka::KA{DIM,NDF}) where{DIM,NDF}
    trees = ka.kdata.field.trees
    kinfo = ka.kinfo

    # Find max vs_num across all local cells that will be processed.
    max_vsn = 1
    for i in eachindex(trees.data)
        for ps_data in trees.data[i]
            isa(ps_data, InsideSolidData) && continue
            ps_data.bound_enc < 0 && continue
            max_vsn = max(max_vsn, ps_data.vs_data.vs_num)
        end
    end

    ws_sL  = zeros(Float64, max_vsn, NDF)
    ws_sR  = zeros(Float64, max_vsn, NDF)
    ws_swL = zeros(Float64, DIM + 2)
    ws_swR = zeros(Float64, DIM + 2)

    @inbounds for i in eachindex(trees.data)
        @inbounds for j in eachindex(trees.data[i])
            ps_data = trees.data[i][j]
            isa(ps_data,InsideSolidData) && continue
            ps_data.bound_enc<0 && continue # solid_cells
            neighbor = ps_data.neighbor
            for dir = 1:DIM
                iL = 2 * dir - 1
                iR = 2 * dir
                update_slope!(
                    Val(neighbor.state[iL]),
                    Val(neighbor.state[iR]),
                    ps_data,
                    kinfo,
                    neighbor.data[iL],
                    neighbor.data[iR],
                    dir,
                    ws_sL, ws_sR, ws_swL, ws_swR,
                )
            end
        end
    end
end

"""
$(TYPEDSIGNATURES)
Physical-space refinement level of `ps_data`, reconstructed from its `ds[1]`.
A tree spans `(geometry[2]-geometry[1])/trees_num[1]` along axis 1, and an
AMR cell at level `L` has `ds[1] = tree_extent / 2^L`.
"""
@inline function cell_level(ps_data::PsData, kinfo::KInfo)
    cfg = kinfo.config
    tree_extent = (cfg.geometry[2] - cfg.geometry[1]) / cfg.trees_num[1]
    Int(round(log2(tree_extent / ps_data.ds[1])))
end

"""
$(TYPEDSIGNATURES)
Per-level slope pass with transverse correction. For every local cell whose
level equals `L`:

  1. Compute `sw`/`sdf` in **every** direction with the original dispatcher
     (same kernels as [`update_slope!`](@ref)) — this initialises the
     non-Val{-1} faces (same-level or hanging-fine sides) which the
     transverse pass alone does not touch.
  2. Re-compute, with the transverse-projection variants, only the faces
     that actually have a coarse (`Val{-1}`) neighbour, where the spurious
     transverse contribution lives.

This must be invoked from the coarsest level upward, with
[`slope_exchange_level!`](@ref) called after every level so the corrected
coarse-side `sw`/`sdf` are visible across ranks before the next level reads
them.
"""
function update_slope_transverse_level!(ka::KA{DIM,NDF}, L::Integer) where{DIM,NDF}
    trees = ka.kdata.field.trees
    kinfo = ka.kinfo

    max_vsn = 1
    for i in eachindex(trees.data)
        for ps_data in trees.data[i]
            isa(ps_data, InsideSolidData) && continue
            ps_data.bound_enc < 0 && continue
            cell_level(ps_data, kinfo) == L || continue
            max_vsn = max(max_vsn, ps_data.vs_data.vs_num)
        end
    end
    ws_sL  = zeros(Float64, max_vsn, NDF)
    ws_sR  = zeros(Float64, max_vsn, NDF)
    ws_swL = zeros(Float64, DIM + 2)
    ws_swR = zeros(Float64, DIM + 2)

    @inbounds for i in eachindex(trees.data)
        @inbounds for j in eachindex(trees.data[i])
            ps_data = trees.data[i][j]
            isa(ps_data, InsideSolidData) && continue
            ps_data.bound_enc < 0 && continue
            cell_level(ps_data, kinfo) == L || continue
            neighbor = ps_data.neighbor

            # (1) Baseline: all directions via the original dispatcher.
            for dir = 1:DIM
                iL = 2 * dir - 1
                iR = 2 * dir
                sL_state = neighbor.state[iL]
                sR_state = neighbor.state[iR]
                (sL_state == -1 || sR_state == -1) && continue
                update_slope!(
                    Val(sL_state),
                    Val(sR_state),
                    ps_data,
                    kinfo,
                    neighbor.data[iL],
                    neighbor.data[iR],
                    dir,
                    ws_sL, ws_sR, ws_swL, ws_swR,
                )
            end

            # (2) Transverse correction on the coarse-neighbour faces. Three
            #     sub-cases per direction:
            #       (a) both sides have neighbour lists and at least one side
            #           is Val{-1}             → two-sided inner-transverse
            #       (b) sR_state == 0 (domain boundary on the R side) and
            #           sL_state == -1         → one-sided bound-transverse using L
            #       (c) sL_state == 0 and sR_state == -1
            #                                  → one-sided bound-transverse using R
            #     Faces with no coarse side (no Val{-1}) keep baseline. Faces
            #     whose non-boundary side is a solid neighbour (bound_enc<0)
            #     also keep baseline — they go through update_slope_bound_*
            #     against the fluid side and the transverse projection on a
            #     solid cell has no well-defined sw to project with.
            DIM == 1 && continue
            for dir = 1:DIM
                iL = 2 * dir - 1
                iR = 2 * dir
                sL_state = neighbor.state[iL]
                sR_state = neighbor.state[iR]
                Ldata = neighbor.data[iL]
                Rdata = neighbor.data[iR]
                if sL_state != 0 && sR_state != 0
                    # case (a): both sides have neighbours.
                    (sL_state == -1 || sR_state == -1) || continue
                    # (Ldata[1].bound_enc < 0 || Rdata[1].bound_enc < 0) && continue
                    dsL = ps_data.midpoint[dir] - Ldata[1].midpoint[dir]
                    dsR = ps_data.midpoint[dir] - Rdata[1].midpoint[dir]
                    update_slope_inner_vs_transverse!(ps_data, Ldata, Rdata, dsL, dsR, dir, ws_sL, ws_sR)
                    if Ldata[1].bound_enc<0
                        update_slope_bound_vs_transverse!(ps_data,Rdata,dsR,dir,ws_sL)
                    elseif Rdata[1].bound_enc<0
                        update_slope_bound_vs_transverse!(ps_data,Ldata,dsL,dir,ws_sL)
                    end
                    # update_slope_inner_ps_transverse!(ps_data, Ldata, Rdata, dsL, dsR, dir, ws_swL, ws_swR)
                elseif sR_state == 0 && sL_state == -1
                    # case (b)
                    Ldata[1].bound_enc < 0 && continue
                    dsL = ps_data.midpoint[dir] - Ldata[1].midpoint[dir]
                    update_slope_bound_vs_transverse!(ps_data, Ldata, dsL, dir, ws_sL)
                    # update_slope_bound_ps_transverse!(ps_data, Ldata, dsL, dir, ws_swL)
                elseif sL_state == 0 && sR_state == -1
                    # case (c)
                    Rdata[1].bound_enc < 0 && continue
                    dsR = ps_data.midpoint[dir] - Rdata[1].midpoint[dir]
                    update_slope_bound_vs_transverse!(ps_data, Rdata, dsR, dir, ws_sL)
                    # update_slope_bound_ps_transverse!(ps_data, Rdata, dsR, dir, ws_swL)
                end
            end
        end
    end
    return nothing
end

"""
$(TYPEDSIGNATURES)
Globally smallest local cell level across all ranks. Used by [`slope!`](@ref)
to decide the starting point of the coarsest-to-finest transverse-correction
sweep.
"""
function min_cell_level(ka::KA)
    trees = ka.kdata.field.trees
    kinfo = ka.kinfo
    local_min = typemax(Int)
    for i in eachindex(trees.data)
        for ps_data in trees.data[i]
            isa(ps_data, InsideSolidData) && continue
            ps_data.bound_enc < 0 && continue
            L = cell_level(ps_data, kinfo)
            L < local_min && (local_min = L)
        end
    end
    local_min == typemax(Int) && (local_min = 0)
    Int(MPI.Allreduce(local_min, MPI.MIN, MPI.COMM_WORLD))
end

"""
$(TYPEDSIGNATURES)
Single-level variant of [`update_slope!`](@ref): re-computes `sw`/`sdf` for
every local cell whose physical-space level equals `L`, using the original
(non-transverse-corrected) inner kernels. Used by [`slope!`](@ref) to seed
the sweep at `L_min` — those cells have no coarser neighbour and therefore
their first-pass slopes are already final.
"""
function update_slope_level!(ka::KA{DIM,NDF}, L::Integer) where{DIM,NDF}
    trees = ka.kdata.field.trees
    kinfo = ka.kinfo

    max_vsn = 1
    for i in eachindex(trees.data)
        for ps_data in trees.data[i]
            isa(ps_data, InsideSolidData) && continue
            ps_data.bound_enc < 0 && continue
            cell_level(ps_data, kinfo) == L || continue
            max_vsn = max(max_vsn, ps_data.vs_data.vs_num)
        end
    end
    ws_sL  = zeros(Float64, max_vsn, NDF)
    ws_sR  = zeros(Float64, max_vsn, NDF)
    ws_swL = zeros(Float64, DIM + 2)
    ws_swR = zeros(Float64, DIM + 2)

    @inbounds for i in eachindex(trees.data)
        @inbounds for j in eachindex(trees.data[i])
            ps_data = trees.data[i][j]
            isa(ps_data, InsideSolidData) && continue
            ps_data.bound_enc < 0 && continue
            cell_level(ps_data, kinfo) == L || continue
            neighbor = ps_data.neighbor
            for dir = 1:DIM
                iL = 2 * dir - 1
                iR = 2 * dir
                update_slope!(
                    Val(neighbor.state[iL]),
                    Val(neighbor.state[iR]),
                    ps_data,
                    kinfo,
                    neighbor.data[iL],
                    neighbor.data[iR],
                    dir,
                    ws_sL, ws_sR, ws_swL, ws_swR,
                )
            end
        end
    end
    return nothing
end


function update_macro_slope!(ka::KA{DIM}) where{DIM}
    trees = ka.kdata.field.trees
    kinfo = ka.kinfo
    @inbounds for i in eachindex(trees.data)
        @inbounds for j in eachindex(trees.data[i])
            ps_data = trees.data[i][j]
            isa(ps_data, InsideSolidData) && continue
            ps_data.bound_enc < 0 && continue
            vs_data = ps_data.vs_data
            for dir = 1:DIM
                ps_data.sw[:,dir] .= @views calc_w0(vs_data.midpoint,vs_data.sdf[:,:,dir],vs_data.weight,kinfo)
            end
        end
    end
end
function slope!(p4est::P_pxest_t,ka::KA{DIM}) where{DIM}
    L_min = min_cell_level(ka)
    L_max = ka.kinfo.config.solver.AMR_PS_MAXLEVEL

    # Coarsest level: no Val{-1} neighbours, so the original kernels already
    # give the final sdf for this level. Publish it only when finer levels
    # need it for transverse projection; a full final exchange below publishes
    # the post-processed sw/sdf used by flux reconstruction.
    update_slope_level!(ka, L_min)
    L_min <= L_max && slope_exchange_level!(p4est, ka, L_min)

    # Finer levels: each level computes its sw/sdf with the transverse
    # correction applied on coarse-neighbour faces, then exchanges so the
    # next level reads the freshly-corrected coarse data. Each exchange is
    # filtered to mirrors/ghosts whose level equals `L`, so the per-level
    # communication cost scales with the size of that level alone.
    for L in (L_min + 1):L_max
        update_slope_transverse_level!(ka, L)
        slope_exchange_level!(p4est, ka, L)
    end
    update_macro_slope!(ka)
    sw_exchange!(p4est, ka)
    return nothing
end
