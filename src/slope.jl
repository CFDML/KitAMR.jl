function vanleer(sL::AbstractArray, sR::AbstractArray)
    SL = abs.(sL)
    SR = abs.(sR)
    @. (sign(sL) + sign(sR)) * SL * SR / (SL + SR + EPS)
end
function vanleer(sL::Real, sR::Real)
    SL = abs(sL)
    SR = abs(sR)
    (sign(sL) + sign(sR)) * SL * SR / (SL + SR + EPS)
end
function diff_L!(vs_data::AbstractVsData{DIM,NDF}, vs_data_n::AbstractVsData{DIM,NDF}, dsL::Float64, sL::AbstractMatrix) where{DIM,NDF}
    index = 1
    flag = 0.0
    level = vs_data.level
    level_n = vs_data_n.level
    df = vs_data.df
    dfn = vs_data_n.df
    for i = 1:vs_data.vs_num
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
function diff_R!(vs_data::AbstractVsData{DIM,NDF}, vs_data_n::AbstractVsData{DIM,NDF}, dsR::Float64, sR::AbstractMatrix) where{DIM,NDF}
    index = 1
    flag = 0.0
    level = vs_data.level
    level_n = vs_data_n.level
    df = vs_data.df
    dfn = vs_data_n.df
    for i = 1:vs_data.vs_num
        if level[i] == level_n[index]
            for j = 1:NDF
                sR[i, j] -= (df[i, j] - dfn[index, j]) / dsR
            end
            index += 1
        elseif level[i] < level_n[index]
            while flag != 1.0
                for j = 1:NDF
                    sR[i, j] -=
                        (df[i, j] - dfn[index, j]) / 2^(DIM * (level_n[index] - level[i])) /
                        dsR
                end
                flag += 1 / 2^(DIM * (level_n[index] - level[i]))
                index += 1
            end
            flag = 0.0
        else
            for j = 1:NDF
                sR[i, j] -= (df[i, j] - dfn[index, j]) / dsR
            end
            flag += 1 / 2^(DIM * (level[i] - level_n[index]))
            if flag == 1.0
                index += 1
                flag = 0.0
            end
        end
    end
end
function update_slope_inner_vs!(
    vs_data::AbstractVsData{DIM,NDF},
    L_data::AbstractVsData{DIM,NDF},
    R_data::AbstractVsData{DIM,NDF},
    dsL::Float64,
    dsR::Float64,
    dir::Int,
) where{DIM,NDF}
    vs_num = vs_data.vs_num
    sL = zeros(vs_num, NDF)
    sR = zeros(vs_num, NDF)
    diff_L!(vs_data, L_data, dsL, sL)
    diff_R!(vs_data, R_data, dsR, sR)
    vs_data.sdf[:, :, dir] .= vanleer(sL, sR)
end
function update_slope_Lbound_vs!(
    vs_data::AbstractVsData{DIM,NDF},
    R_data::AbstractVsData{DIM,NDF},
    dsR::Float64,
    dir::Int,
) where{DIM,NDF}
    vs_num = vs_data.vs_num
    sR = zeros(vs_num, NDF)
    diff_R!(vs_data, R_data, dsR, sR)
    vs_data.sdf[:, :, dir] .= sR
end
function update_slope_Rbound_vs!(
    vs_data::AbstractVsData{DIM,NDF},
    L_data::AbstractVsData{DIM,NDF},
    dsL::Float64,
    dir::Int,
) where{DIM,NDF}
    vs_num = vs_data.vs_num
    sL = zeros(vs_num, NDF)
    diff_L!(vs_data, L_data, dsL, sL)
    vs_data.sdf[:, :, dir] .= sL
end
function update_slope_inner!(
    ::Val{1},
    ::Val{1},
    ps_data::T1,
    global_data::Global_Data{DIM,NDF},
    Ldata::T2,
    Rdata::T2,
    dir::Integer,
) where {T1<:AbstractPsData,T2<:Array,DIM,NDF}
    ds = ps_data.ds[dir]
    swL = (ps_data.w - Ldata[1].w) / (ds)
    swR = (Rdata[1].w - ps_data.w) / (ds)
    sw = vanleer(swL, swR)
    gradmax = global_data.status.gradmax
    for i in eachindex(gradmax)
       @inbounds gradmax[i] = max(gradmax[i], abs(sw[i]))
    end
    ps_data.sw[:, dir] .= sw
    ds = ps_data.ds[dir]
    update_slope_inner_vs!(ps_data.vs_data, Ldata[1].vs_data, Rdata[1].vs_data, ds, ds, dir)
end
function update_slope_inner!(
    ::Val{0},
    ::Val{1},
    ps_data::T1,
    global_data::Global_Data{DIM,NDF},
    Ldata::T2,
    Rdata::T2,
    dir::Integer,
) where {T1<:AbstractPsData,T2<:Array,DIM,NDF}
    sw = (Rdata[1].w - ps_data.w) / ps_data.ds[dir]
    gradmax = global_data.status.gradmax
    for i in eachindex(gradmax)
       @inbounds gradmax[i] = max(gradmax[i], abs(sw[i]))
    end
    ps_data.sw[:, dir] .= sw
    ds = ps_data.ds[dir]
    vs_data = ps_data.vs_data
    update_slope_Lbound_vs!(vs_data, Rdata[1].vs_data, ds, dir)
    faceid = 2*dir-1
    domain = global_data.config.domain[faceid]
    if isdefined(domain,:bc)
        bc = get_bc(domain.bc)
        sL = discrete_maxwell(vs_data.midpoint,bc,global_data)
        sL .= @. (vs_data.df-sL)/ds
        vs_data.sdf[:,:,dir] .= @views vanleer(sL,vs_data.sdf[:,:,dir])
    end
end
function update_slope_inner!(
    ::Val{1},
    ::Val{0},
    ps_data::T1,
    global_data::Global_Data{DIM,NDF},
    Ldata::T2,
    Rdata::T2,
    dir::Integer,
) where {T1<:AbstractPsData,T2<:Array,DIM,NDF}
    sw = (ps_data.w - Ldata[1].w) ./ ps_data.ds[dir]
    gradmax = global_data.status.gradmax
    for i in eachindex(gradmax)
       @inbounds gradmax[i] = max(gradmax[i], abs(sw[i]))
    end
    ps_data.sw[:, dir] .= sw
    ds = ps_data.ds[dir]
    vs_data = ps_data.vs_data
    update_slope_Rbound_vs!(vs_data, Ldata[1].vs_data, ds, dir)
    faceid = 2*dir
    domain = global_data.config.domain[faceid]
    if isdefined(domain,:bc)
        bc = get_bc(domain.bc)
        sR = discrete_maxwell(vs_data.midpoint,bc,global_data)
        sR .= @. (sR-vs_data.df)/ds
        vs_data.sdf[:,:,dir] .= @views vanleer(vs_data.sdf[:,:,dir],sR)
    end
end
function update_slope_Rbound_vs!(
    vs_data::AbstractVsData{DIM,NDF},
    L_datas::Vector,
    dsL::Float64,
    dir::Int,
) where{DIM,NDF}
    sL = zeros(Float64, vs_data.vs_num, NDF)
    for j = 1:2^(DIM-1)
        L_data = L_datas[j].vs_data
        diff_L!(vs_data, L_data, 0.75 * dsL, sL)
    end
    vs_data.sdf[:, :, dir] .= sL / 2^(DIM - 1)
end
function update_slope_Lbound_vs!(
    vs_data::AbstractVsData{DIM,NDF},
    R_datas::Vector,
    dsR::Float64,
    dir::Int,
) where{DIM,NDF}
    sR = zeros(Float64, vs_data.vs_num, NDF)
    for j = 1:2^(DIM-1)
        R_data = R_datas[j].vs_data
        diff_R!(vs_data, R_data, 0.75 * dsR, sR)
    end
    vs_data.sdf[:, :, dir] .= sR / 2^(DIM - 1)
end
function update_slope_inner!(
    ::NeighborNum,
    ::Val{0},
    ps_data::T1,
    global_data::Global_Data{DIM,NDF},
    Ldata::T2,
    ::T2,
    dir::Integer,
) where {T1<:AbstractPsData,T2<:Array,DIM,NDF}
    sw = zeros(Float64, DIM + 2)
    for i = 1:2^(DIM-1)
        sw .+= (ps_data.w - Ldata[i].w)
    end
    sw ./= 2^(DIM - 1) * 0.75 * ps_data.ds[dir]
    gradmax = global_data.status.gradmax
    for i in eachindex(gradmax)
       @inbounds gradmax[i] = max(gradmax[i], abs(sw[i]))
    end
    ps_data.sw[:, dir] .= sw
    ds = ps_data.ds[dir]
    update_slope_Rbound_vs!(ps_data.vs_data, Ldata, ds, dir)
end
function update_slope_inner!(
    ::Val{0},
    ::NeighborNum,
    ps_data::T1,
    global_data::Global_Data{DIM,NDF},
    ::T2,
    Rdata::T2,
    dir::Integer,
) where {T1<:AbstractPsData,T2<:Array,DIM,NDF}
    sw = zeros(Float64, DIM + 2)
    for i = 1:2^(DIM-1)
        sw .+= (Rdata[i].w - ps_data.w)
    end
    sw ./= 2^(DIM - 1) * 0.75 * ps_data.ds[dir]
    gradmax = global_data.status.gradmax
    for i in eachindex(gradmax)
       @inbounds gradmax[i] = max(gradmax[i], abs(sw[i]))
    end
    ps_data.sw[:, dir] .= sw
    ds = ps_data.ds[dir]
    update_slope_Lbound_vs!(ps_data.vs_data, Rdata, ds, dir)
end

function update_slope_inner!(
    ::Val{0},
    ::Val{-1},
    ps_data::T1,
    global_data::Global_Data{DIM,NDF},
    Ldata::T2,
    Rdata::T2,
    dir::Integer,
) where {T1<:AbstractPsData,T2<:Array,DIM,NDF}
    sw = (Rdata[1].w - ps_data.w) / ps_data.ds[dir] / 1.5
    gradmax = global_data.status.gradmax
    for i in eachindex(gradmax)
       @inbounds gradmax[i] = max(gradmax[i], abs(sw[i]))
    end
    ps_data.sw[:, dir] .= sw
    ds = ps_data.ds[dir]
    update_slope_Lbound_vs!(ps_data.vs_data, Rdata[1].vs_data, 1.5 * ds, dir)
end
function update_slope_inner!(
    ::Val{-1},
    ::Val{0},
    ps_data::T1,
    global_data::Global_Data{DIM,NDF},
    Ldata::T2,
    Rdata::T2,
    dir::Integer,
) where {T1<:AbstractPsData,T2<:Array,DIM,NDF}
    sw = (ps_data.w - Ldata[1].w) ./ ps_data.ds[dir] / 1.5
    gradmax = global_data.status.gradmax
    for i in eachindex(gradmax)
       @inbounds gradmax[i] = max(gradmax[i], abs(sw[i]))
    end
    ps_data.sw[:, dir] .= sw
    ds = ps_data.ds[dir]
    update_slope_Rbound_vs!(ps_data.vs_data, Ldata[1].vs_data, 1.5 * ds, dir)
end
function update_slope_inner_vs!(
    vs_data::AbstractVsData{DIM,NDF},
    L_datas::Vector,
    R_datas::Vector,
    dsL::Float64,
    dsR::Float64,
    dir::Int,
) where{DIM,NDF}
    sL = zeros(Float64, vs_data.vs_num, NDF)
    sR = zeros(Float64, vs_data.vs_num, NDF)
    nL = length(L_datas)
    nR = length(R_datas)
    for j = 1:nL
        L_data = L_datas[j].vs_data
        diff_L!(vs_data, L_data, dsL, sL)
    end
    for j = 1:nR
        R_data = R_datas[j].vs_data
        diff_R!(vs_data, R_data, dsR, sR)
    end
    vs_data.sdf[:, :, dir] .= vanleer(sL / nL, sR / nR)
end
function update_slope_inner!(
    ::Val{1},
    ::NeighborNum,
    ps_data::T1,
    global_data::Global_Data{DIM,NDF},
    Ldata::T2,
    Rdata::T2,
    dir::Integer,
) where {T1<:AbstractPsData,T2<:Array,DIM,NDF}
    wR = zeros(DIM + 2)
    for i = 1:2^(DIM-1)
        wR += Rdata[i].w
    end
    ds = ps_data.ds[dir]
    swL = (ps_data.w - Ldata[1].w) / (ds)
    swR = (wR / 2^(DIM - 1) - ps_data.w) / (0.75 * ds)
    sw = vanleer(swL, swR)
    gradmax = global_data.status.gradmax
    for i in eachindex(gradmax)
       @inbounds gradmax[i] = max(gradmax[i], abs(sw[i]))
    end
    ps_data.sw[:, dir] .= sw
    ds = ps_data.ds[dir]
    update_slope_inner_vs!(ps_data.vs_data, Ldata, Rdata, ds, 0.75 * ds, dir)
end
function update_slope_inner!(
    ::NeighborNum,
    ::Val{1},
    ps_data::T1,
    global_data::Global_Data{DIM,NDF},
    Ldata::T2,
    Rdata::T2,
    dir::Integer,
) where {T1<:AbstractPsData,T2<:Array,DIM,NDF}
    wL = zeros(DIM + 2)
    for i = 1:2^(DIM-1)
        wL += Ldata[i].w
    end
    ds = ps_data.ds[dir]
    swL = (ps_data.w - wL / 2^(DIM - 1)) / (0.75 * ds)
    swR = (Rdata[1].w - ps_data.w) / ds
    sw = vanleer(swL, swR)
    gradmax = global_data.status.gradmax
    for i in eachindex(gradmax)
       @inbounds gradmax[i] = max(gradmax[i], abs(sw[i]))
    end
    ps_data.sw[:, dir] .= sw
    ds = ps_data.ds[dir]
    update_slope_inner_vs!(ps_data.vs_data, Ldata, Rdata, 0.75 * ds, ds, dir)
end
function update_slope_inner!(
    ::NeighborNum,
    ::NeighborNum,
    ps_data::T1,
    global_data::Global_Data{DIM,NDF},
    Ldata::T2,
    Rdata::T2,
    dir::Integer,
) where {T1<:AbstractPsData,T2<:Array,DIM,NDF}
    wL = zeros(DIM + 2)
    wR = zeros(DIM + 2)
    for i = 1:2^(DIM-1)
        wL += Ldata[i].w
        wR += Rdata[i].w
    end
    ds = ps_data.ds[dir]
    swL = (ps_data.w - wL / 2^(DIM - 1)) / (0.75 * ds)
    swR = (wR / 2^(DIM - 1) - ps_data.w) / (0.75 * ds)
    sw = vanleer(swL, swR)
    gradmax = global_data.status.gradmax
    for i in eachindex(gradmax)
       @inbounds gradmax[i] = max(gradmax[i], abs(sw[i]))
    end
    ps_data.sw[:, dir] .= sw
    ds = ps_data.ds[dir]
    update_slope_inner_vs!(ps_data.vs_data, Ldata, Rdata, 0.75 * ds, 0.75 * ds, dir)
end
function update_slope_inner!(
    ::Val{1},
    ::Val{-1},
    ps_data::T1,
    global_data::Global_Data{DIM,NDF},
    Ldata::T2,
    Rdata::T2,
    dir::Integer,
) where {T1<:AbstractPsData,T2<:Array,DIM,NDF}
    ds = ps_data.ds[dir]
    swL = (ps_data.w - Ldata[1].w) / ds
    swR = (Rdata[1].w - ps_data.w) / (1.5 * ds)
    sw = vanleer(swL, swR)
    gradmax = global_data.status.gradmax
    for i in eachindex(gradmax)
       @inbounds gradmax[i] = max(gradmax[i], abs(sw[i]))
    end
    ps_data.sw[:, dir] .= sw
    ds = ps_data.ds[dir]
    update_slope_inner_vs!(ps_data.vs_data, Ldata, Rdata, ds, 1.5 * ds, dir)
end
function update_slope_inner!(
    ::Val{-1},
    ::Val{1},
    ps_data::T1,
    global_data::Global_Data{DIM,NDF},
    Ldata::T2,
    Rdata::T2,
    dir::Integer,
) where {T1<:AbstractPsData,T2<:Array,DIM,NDF}
    ds = ps_data.ds[dir]
    swL = (ps_data.w - Ldata[1].w) / (1.5 * ds)
    swR = (Rdata[1].w - ps_data.w) / ds
    sw = vanleer(swL, swR)
    gradmax = global_data.status.gradmax
    for i in eachindex(gradmax)
       @inbounds gradmax[i] = max(gradmax[i], abs(sw[i]))
    end
    ps_data.sw[:, dir] .= sw
    ds = ps_data.ds[dir]
    update_slope_inner_vs!(ps_data.vs_data, Ldata, Rdata, 1.5 * ds, ds, dir)
end
function update_slope_inner!(
    ::Val{-1},
    ::Val{-1},
    ps_data::T1,
    global_data::Global_Data{DIM,NDF},
    Ldata::T2,
    Rdata::T2,
    dir::Integer,
) where {T1<:AbstractPsData,T2<:Array,DIM,NDF}
    ds = ps_data.ds[dir]
    swL = (ps_data.w - Ldata[1].w) / (1.5 * ds)
    swR = (Rdata[1].w - ps_data.w) / (1.5 * ds)
    sw = vanleer(swL, swR)
    gradmax = global_data.status.gradmax
    for i in eachindex(gradmax)
       @inbounds gradmax[i] = max(gradmax[i], abs(sw[i]))
    end
    ps_data.sw[:, dir] .= sw
    ds = ps_data.ds[dir]
    update_slope_inner_vs!(ps_data.vs_data, Ldata, Rdata, 1.5 * ds, 1.5 * ds, dir)
end
function update_slope_inner!(
    ::NeighborNum,
    ::Val{-1},
    ps_data::T1,
    global_data::Global_Data{DIM,NDF},
    Ldata::T2,
    Rdata::T2,
    dir::Integer,
) where {T1<:AbstractPsData,T2<:Array,DIM,NDF}
    wL = zeros(DIM + 2)
    for i = 1:2^(DIM-1)
        wL += Ldata[i].w
    end
    ds = ps_data.ds[dir]
    swL = (ps_data.w - wL / 2^(DIM - 1)) / (0.75 * ds)
    swR = (Rdata[1].w - ps_data.w) / (1.5 * ds)
    sw = vanleer(swL, swR)
    gradmax = global_data.status.gradmax
    for i in eachindex(gradmax)
       @inbounds gradmax[i] = max(gradmax[i], abs(sw[i]))
    end
    ps_data.sw[:, dir] .= sw
    ds = ps_data.ds[dir]
    update_slope_inner_vs!(ps_data.vs_data, Ldata, Rdata, 0.75 * ds, 1.5 * ds, dir)
end
function update_slope_inner!(
    ::Val{-1},
    ::NeighborNum,
    ps_data::T1,
    global_data::Global_Data{DIM,NDF},
    Ldata::T2,
    Rdata::T2,
    dir::Integer,
) where {T1<:AbstractPsData,T2<:Array,DIM,NDF}
    wR = zeros(DIM + 2)
    for i = 1:2^(DIM-1)
        wR += Rdata[i].w
    end
    ds = ps_data.ds[dir]
    swL = (ps_data.w - Ldata[1].w) / (1.5 * ds)
    swR = (wR / 2^(DIM - 1) - ps_data.w) / (0.75 * ds)
    sw = vanleer(swL, swR)
    gradmax = global_data.status.gradmax
    for i in eachindex(gradmax)
       @inbounds gradmax[i] = max(gradmax[i], abs(sw[i]))
    end
    ps_data.sw[:, dir] .= sw
    ds = ps_data.ds[dir]
    update_slope_inner_vs!(ps_data.vs_data, Ldata, Rdata, 1.5 * ds, 0.75 * ds, dir)
end
function update_slope!(amr::AMR{DIM,NDF}) where{DIM,NDF}
    trees = amr.field.trees
    global_data = amr.global_data
    @inbounds for i in eachindex(trees.data)
        @inbounds for j in eachindex(trees.data[i])
            ps_data = trees.data[i][j]
            isa(ps_data,InsideSolidData) && continue
            ps_data.bound_enc<0 && continue # solid_cells
            neighbor = ps_data.neighbor
            for dir = 1:DIM
                iL = 2 * dir - 1
                iR = 2 * dir
                update_slope_inner!(
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
