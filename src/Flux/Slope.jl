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
function minmod(sL::Real, sR::Real)
    SL = abs(sL)
    SR = abs(sR)
    0.5(sign(sL)+sign(sR))*min(SL, SR)
end
"""
$(TYPEDSIGNATURES)
Difference through all velocity cells in two neighboring physical cells.
"""
function diff_vs!(
    vs_data::AbstractVsData{DIM,NDF},
    vs_data_n::AbstractVsData{DIM,NDF},
    dsL::Float64,
    sL::AbstractMatrix,
) where {DIM,NDF}
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
) where {DIM,NDF}
    sdf = zeros(Float64, vs_data.vs_num, NDF)
    N = length(neighbor_datas)
    for j = 1:N
        neighbor_data = neighbor_datas[j].vs_data
        diff_vs!(vs_data, neighbor_data, ds, sdf)
    end
    vs_data.sdf[:, :, dir] .= sdf / N
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
) where {DIM,NDF}
    sL = zeros(Float64, vs_data.vs_num, NDF)
    sR = zeros(Float64, vs_data.vs_num, NDF)
    nL = length(Ldata);
    nR = length(Rdata)
    for j = 1:nL
        L_data = Ldata[j].vs_data
        diff_vs!(vs_data, L_data, dsL, sL)
    end
    for j = 1:nR
        R_data = Rdata[j].vs_data
        diff_vs!(vs_data, R_data, dsR, sR)
    end
    vs_data.sdf[:, :, dir] .= vanleer(sL / nL, sR / nR)
    return nothing
end
"""
$(TYPEDSIGNATURES)
"""
function update_slope_inner_ps!(
    ps_data::PS_Data{DIM,NDF},
    Ldata::AbstractVector,
    Rdata::AbstractVector,
    dsL::Float64,
    dsR::Float64,
    dir::Int,
) where {DIM,NDF}
    swL = zeros(Float64, DIM+2);
    swR = zeros(Float64, DIM+2)
    nL = length(Ldata);
    nR = length(Rdata)
    for j = 1:nL
        @. swL += (ps_data.w-Ldata[j].w)/dsL
    end
    for j = 1:nR
        @. swR += (ps_data.w-Rdata[j].w)/dsR
    end
    swL/=nL;
    swR/=nR
    for j in eachindex(swL)
        ps_data.sw[j, dir] = abs(swL[j])>abs(swR[j]) ? swL[j] : swR[j]
    end
    return nothing
end
"""
$(TYPEDSIGNATURES)
"""
function update_slope_bound_ps!(
    ps_data::PS_Data{DIM,NDF},
    Ldata::AbstractVector,
    dsL::Float64,
    dir::Int,
) where {DIM,NDF}
    swL = zeros(Float64, DIM+2)
    nL = length(Ldata)
    for j = 1:nL
        @. swL += (ps_data.w-Ldata[j].w)/dsL
    end
    swL/=nL
    ps_data.sw[:, dir] .= swL
    return nothing
end


# Inner
"""
$(TYPEDSIGNATURES)
"""
function update_slope!(
    ::Val{1},
    ::Val{1},
    ps_data::PS_Data,
    global_data::Global_Data{DIM,NDF},
    Ldata::AbstractVector,
    Rdata::AbstractVector,
    dir::Integer,
) where {DIM,NDF}
    if Ldata[1].bound_enc<0&&Rdata[1].bound_enc<0
        ps_data.vs_data.sdf[:, :, dir] .= 0.0
    elseif Ldata[1].bound_enc<0
        ds = ps_data.midpoint[dir]-Rdata[1].midpoint[dir]
        update_slope_bound_vs!(ps_data.vs_data, Rdata, ds, dir)
        update_slope_bound_ps!(ps_data, Rdata, ds, dir)
    elseif Rdata[1].bound_enc<0
        ds = ps_data.midpoint[dir]-Ldata[1].midpoint[dir]
        update_slope_bound_vs!(ps_data.vs_data, Ldata, ds, dir)
        update_slope_bound_ps!(ps_data, Ldata, ds, dir)
    else
        ds = ps_data.ds[dir]
        update_slope_inner_vs!(ps_data.vs_data, Ldata, Rdata, ds, -ds, dir)
        update_slope_inner_ps!(ps_data, Ldata, Rdata, ds, -ds, dir)
    end
    return nothing
end
"""
$(TYPEDSIGNATURES)
"""
function update_slope!(
    ::Val{1},
    ::NeighborNum,
    ps_data::PS_Data,
    global_data::Global_Data{DIM,NDF},
    Ldata::AbstractVector,
    Rdata::AbstractVector,
    dir::Integer,
) where {DIM,NDF}
    ds = ps_data.ds[dir]
    update_slope_inner_vs!(ps_data.vs_data, Ldata, Rdata, ds, -0.75 * ds, dir)
    update_slope_inner_ps!(ps_data, Ldata, Rdata, ds, -0.75 * ds, dir)
end
"""
$(TYPEDSIGNATURES)
"""
function update_slope!(
    ::NeighborNum,
    ::Val{1},
    ps_data::PS_Data,
    global_data::Global_Data{DIM,NDF},
    Ldata::AbstractVector,
    Rdata::AbstractVector,
    dir::Integer,
) where {DIM,NDF}
    ds = ps_data.ds[dir]
    update_slope_inner_vs!(ps_data.vs_data, Ldata, Rdata, 0.75 * ds, -ds, dir)
    update_slope_inner_ps!(ps_data, Ldata, Rdata, 0.75 * ds, -ds, dir)
end
"""
$(TYPEDSIGNATURES)
"""
function update_slope!(
    ::NeighborNum,
    ::NeighborNum,
    ps_data::PS_Data,
    global_data::Global_Data{DIM,NDF},
    Ldata::AbstractVector,
    Rdata::AbstractVector,
    dir::Integer,
) where {DIM,NDF}
    ds = ps_data.ds[dir]
    update_slope_inner_vs!(ps_data.vs_data, Ldata, Rdata, 0.75 * ds, -0.75 * ds, dir)
    update_slope_inner_ps!(ps_data, Ldata, Rdata, 0.75 * ds, -0.75 * ds, dir)
end
"""
$(TYPEDSIGNATURES)
"""
function update_slope!(
    ::Val{-1},
    ::Val{-1},
    ps_data::PS_Data,
    global_data::Global_Data{DIM,NDF},
    Ldata::AbstractVector,
    Rdata::AbstractVector,
    dir::Integer,
) where {DIM,NDF}
    ds = ps_data.ds[dir]
    update_slope_inner_vs!(ps_data.vs_data, Ldata, Rdata, 1.5 * ds, -1.5 * ds, dir)
    update_slope_inner_ps!(ps_data, Ldata, Rdata, 1.5 * ds, -1.5 * ds, dir)
end
"""
$(TYPEDSIGNATURES)
"""
function update_slope!(
    ::Val{1},
    ::Val{-1},
    ps_data::PS_Data,
    global_data::Global_Data{DIM,NDF},
    Ldata::AbstractVector,
    Rdata::AbstractVector,
    dir::Integer,
) where {DIM,NDF}
    ds = ps_data.ds[dir]
    update_slope_inner_vs!(ps_data.vs_data, Ldata, Rdata, ds, -1.5 * ds, dir)
    update_slope_inner_ps!(ps_data, Ldata, Rdata, ds, -1.5 * ds, dir)
end
"""
$(TYPEDSIGNATURES)
"""
function update_slope!(
    ::Val{-1},
    ::Val{1},
    ps_data::PS_Data,
    global_data::Global_Data{DIM,NDF},
    Ldata::AbstractVector,
    Rdata::AbstractVector,
    dir::Integer,
) where {DIM,NDF}
    ds = ps_data.ds[dir]
    update_slope_inner_vs!(ps_data.vs_data, Ldata, Rdata, 1.5 * ds, -ds, dir)
    update_slope_inner_ps!(ps_data, Ldata, Rdata, 1.5 * ds, -ds, dir)
end
"""
$(TYPEDSIGNATURES)
"""
function update_slope!(
    ::NeighborNum,
    ::Val{-1},
    ps_data::PS_Data,
    global_data::Global_Data{DIM,NDF},
    Ldata::AbstractVector,
    Rdata::AbstractVector,
    dir::Integer,
) where {DIM,NDF}
    ds = ps_data.ds[dir]
    update_slope_inner_vs!(ps_data.vs_data, Ldata, Rdata, 0.75 * ds, -1.5 * ds, dir)
    update_slope_inner_ps!(ps_data, Ldata, Rdata, 0.75 * ds, -1.5 * ds, dir)
end
"""
$(TYPEDSIGNATURES)
"""
function update_slope!(
    ::Val{-1},
    ::NeighborNum,
    ps_data::PS_Data,
    global_data::Global_Data{DIM,NDF},
    Ldata::AbstractVector,
    Rdata::AbstractVector,
    dir::Integer,
) where {DIM,NDF}
    ds = ps_data.ds[dir]
    update_slope_inner_vs!(ps_data.vs_data, Ldata, Rdata, 1.5 * ds, -0.75 * ds, dir)
    update_slope_inner_ps!(ps_data, Ldata, Rdata, 1.5 * ds, -0.75 * ds, dir)
end

# Boundary
"""
$(TYPEDSIGNATURES)
"""
function update_slope!(
    ::Val{0},
    ::Val{1},
    ps_data::PS_Data,
    global_data::Global_Data{DIM,NDF},
    Ldata::AbstractVector,
    Rdata::AbstractVector,
    dir::Integer,
) where {DIM,NDF}
    ds = ps_data.midpoint[dir]-Rdata[1].midpoint[dir]
    vs_data = ps_data.vs_data
    update_slope_bound_vs!(vs_data, Rdata, ds, dir)
    update_slope_bound_ps!(ps_data, Rdata, ds, dir)
end
"""
$(TYPEDSIGNATURES)
"""
function update_slope!(
    ::Val{1},
    ::Val{0},
    ps_data::PS_Data,
    global_data::Global_Data{DIM,NDF},
    Ldata::AbstractVector,
    Rdata::AbstractVector,
    dir::Integer,
) where {DIM,NDF}
    ds = ps_data.midpoint[dir]-Ldata[1].midpoint[dir]
    vs_data = ps_data.vs_data
    update_slope_bound_vs!(vs_data, Ldata, ds, dir)
    update_slope_bound_ps!(ps_data, Ldata, ds, dir)
end
"""
$(TYPEDSIGNATURES)
"""
function update_slope!(
    ::NeighborNum,
    ::Val{0},
    ps_data::PS_Data,
    global_data::Global_Data{DIM,NDF},
    Ldata::AbstractVector,
    Rdata::AbstractVector,
    dir::Integer,
) where {DIM,NDF}
    ds = ps_data.midpoint[dir]-Ldata[1].midpoint[dir]
    update_slope_bound_vs!(ps_data.vs_data, Ldata, ds, dir)
    update_slope_bound_ps!(ps_data, Ldata, ds, dir)
end
"""
$(TYPEDSIGNATURES)
"""
function update_slope!(
    ::Val{0},
    ::NeighborNum,
    ps_data::PS_Data,
    global_data::Global_Data{DIM,NDF},
    Ldata::AbstractVector,
    Rdata::AbstractVector,
    dir::Integer,
) where {DIM,NDF}
    ds = ps_data.midpoint[dir]-Rdata[1].midpoint[dir]
    update_slope_bound_vs!(ps_data.vs_data, Rdata, ds, dir)
    update_slope_bound_ps!(ps_data, Rdata, ds, dir)
end
"""
$(TYPEDSIGNATURES)
"""
function update_slope!(
    ::Val{0},
    ::Val{-1},
    ps_data::PS_Data,
    global_data::Global_Data{DIM,NDF},
    Ldata::AbstractVector,
    Rdata::AbstractVector,
    dir::Integer,
) where {DIM,NDF}
    ds = ps_data.midpoint[dir]-Rdata[1].midpoint[dir]
    update_slope_bound_vs!(ps_data.vs_data, Rdata, ds, dir)
    update_slope_bound_ps!(ps_data, Rdata, ds, dir)
end
"""
$(TYPEDSIGNATURES)
"""
function update_slope!(
    ::Val{-1},
    ::Val{0},
    ps_data::PS_Data,
    global_data::Global_Data{DIM,NDF},
    Ldata::AbstractVector,
    Rdata::AbstractVector,
    dir::Integer,
) where {DIM,NDF}
    ds = ps_data.midpoint[dir]-Ldata[1].midpoint[dir]
    update_slope_bound_vs!(ps_data.vs_data, Ldata, ds, dir)
    update_slope_bound_ps!(ps_data, Ldata, ds, dir)
end

"""
$(TYPEDSIGNATURES)
"""
function update_slope!(amr::KitAMR_Data{DIM,NDF}) where {DIM,NDF}
    trees = amr.field.trees
    global_data = amr.global_data
    @inbounds for i in eachindex(trees.data)
        @inbounds for j in eachindex(trees.data[i])
            ps_data = trees.data[i][j]
            isa(ps_data, InsideSolidData) && continue
            ps_data.bound_enc<0 && continue # solid_cells
            neighbor = ps_data.neighbor
            for dir = 1:DIM
                iL = 2 * dir - 1
                iR = 2 * dir
                update_slope!(
                    Val(neighbor.state[iL]),
                    Val(neighbor.state[iR]),
                    ps_data,
                    global_data,
                    neighbor.data[iL],
                    neighbor.data[iR],
                    dir,
                )
            end
        end
    end
end
