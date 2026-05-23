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
        ps_data.sw[j, dir] = vanleer(ws_swL[j],ws_swR[j])
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
        update_slope_bound_ps!(ps_data, Rdata, ds, dir, ws_swL)
    elseif Rdata[1].bound_enc<0
        ds = ps_data.midpoint[dir]-Ldata[1].midpoint[dir]
        update_slope_bound_vs!(ps_data.vs_data, Ldata, ds, dir, ws_sL)
        update_slope_bound_ps!(ps_data, Ldata, ds, dir, ws_swL)
    else
        ds = ps_data.ds[dir]
        update_slope_inner_vs!(ps_data.vs_data, Ldata, Rdata, ds, -ds, dir, ws_sL, ws_sR)
        update_slope_inner_ps!(ps_data, Ldata, Rdata, ds, -ds, dir, ws_swL, ws_swR)
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
    update_slope_inner_ps!(ps_data, Ldata, Rdata, ds, -0.75 * ds, dir, ws_swL, ws_swR)
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
    update_slope_inner_ps!(ps_data, Ldata, Rdata, 0.75 * ds, -ds, dir, ws_swL, ws_swR)
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
    update_slope_inner_ps!(ps_data, Ldata, Rdata, 0.75 * ds, -0.75 * ds, dir, ws_swL, ws_swR)
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
    update_slope_inner_ps!(ps_data, Ldata, Rdata, 1.5 * ds, -1.5 * ds, dir, ws_swL, ws_swR)
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
    update_slope_inner_ps!(ps_data, Ldata, Rdata, ds, -1.5 * ds, dir, ws_swL, ws_swR)
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
    update_slope_inner_ps!(ps_data, Ldata, Rdata, 1.5 * ds, -ds, dir, ws_swL, ws_swR)
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
    update_slope_inner_ps!(ps_data, Ldata, Rdata, 0.75 * ds, -1.5 * ds, dir, ws_swL, ws_swR)
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
    update_slope_inner_ps!(ps_data, Ldata, Rdata, 1.5 * ds, -0.75 * ds, dir, ws_swL, ws_swR)
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
    update_slope_bound_ps!(ps_data, Rdata, ds, dir, ws_swL)
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
    update_slope_bound_ps!(ps_data, Ldata, ds, dir, ws_swL)
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
    update_slope_bound_ps!(ps_data, Ldata, ds, dir, ws_swL)
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
    update_slope_bound_ps!(ps_data, Rdata, ds, dir, ws_swL)
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
    update_slope_bound_ps!(ps_data, Rdata, ds, dir, ws_swL)
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
    update_slope_bound_ps!(ps_data, Ldata, ds, dir, ws_swL)
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
function slope!(p4est::P_pxest_t,ka::KA)
    update_slope!(ka)
    slope_exchange!(p4est, ka)
    return nothing
end