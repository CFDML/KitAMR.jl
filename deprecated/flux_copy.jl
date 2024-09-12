function calc_face_area(ps_data::AbstractPsData{2}, DIR::Integer)
    return ps_data.ds[FAT[1][DIR]]
end
function calc_face_area(ps_data::AbstractPsData{3}, DIR::Integer)
    return reduce(*, @view(ps_data.ds[FAT[2][DIR]]))
end
function calc_flux!(::BoundaryNeighbor, face::Face{BoundaryFace}, amr::AMR{2,2})
    ps_data = face.data
    global_data = amr.global_data
    faceid = face.faceid
    DIR = get_dir(faceid)
    ROT = get_rot(faceid)
    vs_data = ps_data.vs_data
    midpoint = vs_data.midpoint
    Θ = [heaviside(ROT * midpoint[i, DIR]) for i in axes(midpoint, 1)]
    vn = @views midpoint[:, DIR]
    prim0 = global_data.config.bc[:, faceid]
    df0 = @. vs_data.df - 0.5 * ROT * ps_data.ds[DIR] * @view(vs_data.sdf[:, :, DIR])
    prim0[1] = calc_ρw(vs_data, df0, prim0, Θ, vn)
    F0 = discrete_maxwell(vs_data.midpoint, prim0, global_data)
    F = @. F0 * Θ + df0 * (1.0 - Θ)
    fw = calc_fwb(vs_data, F, vn)
    dsf = calc_face_area(ps_data, DIR)
    Δt = global_data.status.Δt
    @. ps_data.flux += ROT * fw * ds * Δt
    @. vs_data.flux[:, 1] += vn * ROT * dsf * Δt * @view(F[:, 1])
    @. vs_data.flux[:, 2] += vn * ROT * dsf * Δt * @view(F[:, 2])
end
function calc_flux!(::BoundaryNeighbor, face::Face{BoundaryFace}, amr::AMR{3,1})
    ps_data = face.data
    global_data = amr.global_data
    faceid = face.faceid
    DIR = get_dir(faceid)
    ROT = get_rot(faceid)
    vs_data = ps_data.vs_data
    midpoint = vs_data.midpoint
    Θ = [heaviside(ROT * midpoint[i, DIR]) for i in axes(midpoint, 1)]
    vn = @views midpoint[:, DIR]
    prim0 = global_data.config.bc[:, faceid]
    df0 = @. vs_data.df - 0.5 * ROT * ps_data.ds[DIR] * @view(vs_data.sdf[:, :, DIR])
    prim0[1] = calc_ρw(vs_data, df0, prim0, Θ, vn)
    F0 = discrete_maxwell(vs_data.midpoint, prim0, global_data)
    F = @. F0 * Θ + df0 * (1.0 - Θ)
    fw = calc_fwb(vs_data, F, vn)
    dsf = calc_face_area(ps_data, DIR)
    Δt = global_data.status.Δt
    @. ps_data.flux += ROT * fw * dsf * Δt
    @. vs_data.flux += vn * ROT * dsf * Δt * F
end
function update_vs_flux!(
    micro_flux::AbstractMatrix,
    bit_L::AbstractVector,
    vs_data::AbstractVsData{DIM},
    vs_data_n::VS_Data,
    offset::Int,
    ROT::Float64,
) where {DIM}
    index = j = 1
    index_n = offset + 1
    flag = 0.0
    level = vs_data.level
    level_n = vs_data_n.level
    flux = vs_data.flux
    flux_n = vs_data_n.flux
    for i = 1:vs_data.vs_num
        if bit_L[i]
            @. flux[i, :] += @views ROT * micro_flux[index, :]
            if level[i] == level_n[j]
                @. flux_n[j, :] -= @views ROT * micro_flux[index, :]
                j += 1
            elseif level[i] < level_n[j]
                while flag != 1.0
                    @. flux_n[j, :] -= @views ROT * micro_flux[index, :]
                    flag += 1 / 2^(DIM * (level_n[j] - level[i]))
                    j += 1
                end
                flag = 0.0
            else
                @. flux_n[j, :] -=
                    ROT * @view(micro_flux[index, :]) / 2^(DIM * (level[i] - level_n[j]))
                flag += 1 / 2^(DIM * (level[i] - level_n[j]))
                if flag == 1.0
                    j += 1
                    flag = 0.0
                end
            end
            index += 1
        else
            if level[i] == level_n[j]
                @. flux_n[j, :] -= @views ROT * micro_flux[index_n, :]
                @. flux[i, :] += @views ROT * micro_flux[index_n, :]
                j += 1
                index_n += 1
            elseif level[i] < level_n[j]
                while flag != 1.0
                    @. flux_n[j, :] -= @views ROT * micro_flux[index_n, :]
                    @. flux[i, :] +=
                        ROT * @view(micro_flux[index_n, :]) /
                        2^(DIM * (level_n[j] - level[i]))
                    flag += 1 / 2^(DIM * (level_n[j] - level[i]))
                    j += 1
                    index_n += 1
                end
                flag = 0.0
            else
                @. flux[i, :] += @views ROT * micro_flux[index_n, :]
                flag += 1 / 2^(DIM * (level[i] - level_n[j]))
                if flag == 1.0
                    @. flux_n[j, :] -= @views ROT * micro_flux[index_n, :]
                    j += 1
                    index_n += 1
                    flag = 0.0
                end
            end
        end
    end
end
function update_vs_flux!(
    micro_flux::AbstractMatrix,
    bit_L::AbstractVector,
    vs_data::AbstractVsData{DIM},
    vs_data_n::Ghost_VS_Data,
    offset::Int,
    ROT::Float64,
) where {DIM}
    index = j = 1
    index_n = offset + 1
    flag = 0.0
    level = vs_data.level
    level_n = vs_data_n.level
    flux = vs_data.flux
    for i = 1:vs_data.vs_num
        if bit_L[i]
            @. flux[i, :] += @views ROT * micro_flux[index, :]
            if level[i] == level_n[j]
                j += 1
            elseif level[i] < level_n[j]
                while flag != 1.0
                    flag += 1 / 2^(DIM * (level_n[j] - level[i]))
                    j += 1
                end
                flag = 0.0
            else
                flag += 1 / 2^(DIM * (level[i] - level_n[j]))
                if flag == 1.0
                    j += 1
                    flag = 0.0
                end
            end
            index += 1
        else
            if level[i] == level_n[j]
                @. flux[i, :] += @views ROT * micro_flux[index_n, :]
                j += 1
                index_n += 1
            elseif level[i] < level_n[j]
                while flag != 1.0
                    @. flux[i, :] +=
                        ROT * @view(micro_flux[index_n, :]) /
                        2^(DIM * (level_n[j] - level[i]))
                    flag += 1 / 2^(DIM * (level_n[j] - level[i]))
                    j += 1
                    index_n += 1
                end
                flag = 0.0
            else
                @. flux[i, :] += ROT * @view(micro_flux[index_n, :])
                flag += 1 / 2^(DIM * (level[i] - level_n[j]))
                if flag == 1.0
                    j += 1
                    index_n += 1
                    flag = 0.0
                end
            end
        end
    end
end

function update_nflux!(nps_data::PS_Data, fw::AbstractVector)
    nps_data.flux .-= fw
end
function update_nflux!(::Ghost_PS_Data, ::AbstractVector)
    return nothing
end

function make_face_data(
    ps_data::PS_Data{2},
    nps_data::AbstractPsData,
    faceid::Int,
    ROT::Float64,
    DIR::Int,
)
    vs_data = ps_data.vs_data
    midpoint_L = vs_data.midpoint
    vs_data_n = nps_data.vs_data
    midpoint_R = vs_data_n.midpoint
    bit_L = [x <= 0.0 for x in ROT .* @view(midpoint_L[:, DIR])]
    bit_R = [x > 0.0 for x in ROT .* @view(midpoint_R[:, DIR])]
    midpoint_n = Vcat(@view(midpoint_L[bit_L, :]), @view(midpoint_R[bit_R, :]))
    if ps_data.neighbor.state[faceid] == 1
        df_L = @view(vs_data.df[bit_L, :])
        df_R = @view(vs_data_n.df[bit_R, :])
    elseif ps_data.neighbor.state[faceid] > 1
        DIR_t = DIR % 2 + 1
        df_L = @. @view(vs_data.df[bit_L, :]) +
           (nps_data.midpoint[DIR_t] - ps_data.midpoint[DIR_t]) *
           @view(vs_data.sdf[bit_L, :, DIR_t])
        df_R = @view(vs_data_n.df[bit_R, :])
    else
        DIR_t = DIR % 2 + 1
        df_L = @view(vs_data.df[bit_L, :])
        df_R = @. @view(vs_data_n.df[bit_R, :]) +
           (ps_data.midpoint[DIR_t] - nps_data.midpoint[DIR_t]) *
           @view(vs_data_n.sdf[bit_R, :, DIR_t])
    end
    df_n = Vcat(df_L, df_R)
    sdf_n = Vcat(@view(vs_data.sdf[bit_L, :, DIR]), @view(vs_data_n.sdf[bit_R, :, DIR]))
    weight_L = @view(vs_data.weight[bit_L])
    offset = length(weight_L)
    weight_n = Vcat(weight_L, @view(vs_data_n.weight[bit_R]))
    (Face_VS_Data{2,NDF}(weight_n, midpoint_n, df_n, sdf_n), offset, bit_L, bit_R)
end

function make_face_data(
    ps_data::PS_Data{3,NDF},
    nps_data::AbstractPsData,
    faceid::Int,
    ROT::Float64,
    DIR::Int,
) where {NDF}
    vs_data = ps_data.vs_data
    midpoint_L = vs_data.midpoint
    vs_data_n = nps_data.vs_data
    midpoint_R = vs_data_n.midpoint
    bit_L = [x <= 0.0 for x in ROT .* @view(midpoint_L[:, DIR])]
    bit_R = [x > 0.0 for x in ROT .* @view(midpoint_R[:, DIR])]
    midpoint_n = Vcat(@view(midpoint_L[bit_L, :]), @view(midpoint_R[bit_R, :]))
    vn = @view(midpoint_n[:, DIR])
    if ps_data.neighbor.state[faceid] == 1
        df_L = @view(vs_data.df[bit_L, :])
        df_R = @view(vs_data_n.df[bit_R, :])
    elseif ps_data.neighbor.state[faceid] > 1
        df_L = @. @view(vs_data.df[bit_L, :]) +
           (nps_data.midpoint[FAT[2][DIR][1]] - ps_data.midpoint[FAT[2][DIR][1]]) *
           @view(vs_data.sdf[bit_L, :, FAT[2][DIR][1]]) +
           (nps_data.midpoint[FAT[2][DIR][2]] - ps_data.midpoint[FAT[2][DIR][2]]) *
           @view(vs_data.sdf[bit_L, :, FAT[2][DIR][2]])
        df_R = @view(vs_data_n.df[bit_R, :])
    else
        df_L = @view(vs_data.df[bit_L, :])
        df_R = @. @view(vs_data_n.df[bit_R, :]) +
           (ps_data.midpoint[FAT[2][DIR][1]] - nps_data.midpoint[FAT[2][DIR][1]]) *
           @view(vs_data_n.sdf[bit_R, :, FAT[2][DIR][1]]) +
           (ps_data.midpoint[FAT[2][DIR][2]] - nps_data.midpoint[FAT[2][DIR][2]]) *
           @view(vs_data_n.sdf[bit_R, :, FAT[2][DIR][2]])
    end
    df_n = Vcat(df_L, df_R)
    sdf_n = Vcat(@view(vs_data.sdf[bit_L, :, DIR]), @view(vs_data_n.sdf[bit_R, :, DIR]))
    weight_L = @view(vs_data.weight[bit_L])
    offset = length(weight_L)
    weight_n = Vcat(weight_L, @view(vs_data_n.weight[bit_R]))
    (Face_VS_Data{3,NDF}(weight_n, midpoint_n, vn, df_n, sdf_n), offset, bit_L, bit_R)
end
function update_face_data_R!(
    ps_data::PS_Data,
    nps_data::AbstractPsData,
    f_vs_data::Face_VS_Data{2},
    offset::Integer,
    faceid::Integer,
    ROT::Real,
    DIR::Integer,
)
    vs_data_n = nps_data.vs_data
    midpoint_R = vs_data_n.midpoint
    bit_R = [x > 0.0 for x in ROT .* @view(midpoint_R[:, DIR])]
    if ps_data.neighbor.state[faceid] < 1
        df_R = @. @view(vs_data_n.df[bit_R, :]) +
           (ps_data.midpoint[DIR_t] - nps_data.midpoint[DIR_t]) *
           @view(vs_data_n.sdf[bit_R, :, DIR_t])
    else
        DIR_t = DIR % 2 + 1
        df_R = @view(vs_data_n.df[bit_R, :])
    end
    f_vs_data.midpoint =
        Vcat(@view(f_vs_data.midpoint[1:offset, :]), @view(midpoint_R[bit_R, :]))
    f_vs_data.vn = @view(f_vs_data.midpoint[:, DIR])
    f_vs_data.df = Vcat(@view(f_vs_data.df[1:offset, :]), df_R)
    f_vs_data.sdf =
        Vcat(@view(f_vs_data.sdf[1:offset, :]), @view(vs_data_n.sdf[bit_R, :, DIR]))
    f_vs_data.weight =
        Vcat(@view(f_vs_data.weight[1:offset]), @view(vs_data_n.weight[bit_R]))
end
function update_face_data_R!(
    ps_data::PS_Data,
    nps_data::AbstractPsData,
    f_vs_data::Face_VS_Data{3},
    offset::Integer,
    faceid::Integer,
    ROT::Real,
    DIR::Integer,
)
    vs_data_n = nps_data.vs_data
    midpoint_R = vs_data_n.midpoint
    bit_R = [x > 0.0 for x in ROT .* @view(midpoint_R[:, DIR])]
    if ps_data.neighbor.state[faceid] < 1
        df_R = @. @view(vs_data_n.df[bit_R, :]) +
           (ps_data.midpoint[FAT[2][DIR][1]] - nps_data.midpoint[FAT[2][DIR][1]]) *
           @view(vs_data_n.sdf[bit_R, :, FAT[2][DIR][1]]) +
           (ps_data.midpoint[FAT[2][DIR][2]] - nps_data.midpoint[FAT[2][DIR][2]]) *
           @view(vs_data_n.sdf[bit_R, :, FAT[2][DIR][2]])
    else
        df_R = @view(vs_data_n.df[bit_R, :])
    end
    f_vs_data.midpoint =
        Vcat(@view(f_vs_data.midpoint[1:offset, :]), @view(midpoint_R[bit_R, :]))
    f_vs_data.vn = @view(f_vs_data.midpoint[:, DIR])
    f_vs_data.df = Vcat(@view(f_vs_data.df[1:offset, :]), df_R)
    f_vs_data.sdf =
        Vcat(@view(f_vs_data.sdf[1:offset, :]), @view(vs_data_n.sdf[bit_R, :, DIR]))
    f_vs_data.weight =
        Vcat(@view(f_vs_data.weight[1:offset]), @view(vs_data_n.weight[bit_R]))
end
function update_face_data_L!(
    ps_data::PS_Data,
    nps_data::AbstractPsData,
    f_vs_data::Face_VS_Data{2},
    offset::Integer,
    faceid::Integer,
    ROT::Real,
    DIR::Integer,
)
    vs_data = ps_data.vs_data
    midpoint_L = vs_data.midpoint
    bit_L = [x <= 0.0 for x in ROT .* @view(midpoint_L[:, DIR])]
    if ps_data.neighbor.state[faceid] > 1
        DIR_t = DIR % 2 + 1
        df_L = @. @view(vs_data.df[bit_L, :]) +
           (nps_data.midpoint[DIR_t] - ps_data.midpoint[DIR_t]) *
           @view(vs_data.sdf[bit_L, :, DIR_t])
    else
        df_L = @view(vs_data.df[bit_L, :])
    end
    f_vs_data.midpoint =
        Vcat(@view(midpoint_L[bit_L, :]), @view(f_vs_data.midpoint[offset+1:end, :]))
    f_vs_data.vn = @view(f_vs_data.midpoint[:, DIR])
    f_vs_data.df = Vcat(df_L, @view(f_vs_data.df[offset+1:end, :]))
    f_vs_data.sdf =
        Vcat(@view(vs_data.sdf[bit_L, :, DIR]), @view(f_vs_data.sdf[offset+1:end, :, DIR]))
    f_vs_data.weight =
        Vcat(@view(vs_data.weight[bit_L]), @view(f_vs_data.weight[offset+1:end]))
    offset = length(df_L)
    return (bit_L, offset)
end
function update_face_data_L!(
    ps_data::PS_Data,
    nps_data::AbstractPsData,
    f_vs_data::Face_VS_Data{3},
    offset::Integer,
    faceid::Integer,
    ROT::Real,
    DIR::Integer,
)
    vs_data = ps_data.vs_data
    midpoint_L = vs_data.midpoint
    bit_L = [x <= 0.0 for x in ROT .* @view(midpoint_L[:, DIR])]
    if ps_data.neighbor.state[faceid] > 1
        df_L = @. @view(vs_data.df[bit_L, :]) +
           (nps_data.midpoint[FAT[2][DIR][1]] - ps_data.midpoint[FAT[2][DIR][1]]) *
           @view(vs_data.sdf[bit_L, :, FAT[2][DIR][1]]) +
           (nps_data.midpoint[FAT[2][DIR][2]] - ps_data.midpoint[FAT[2][DIR][2]]) *
           @view(vs_data.sdf[bit_L, :, FAT[2][DIR][2]])
    else
        df_L = @view(vs_data.df[bit_L, :])
    end
    f_vs_data.midpoint =
        Vcat(@view(midpoint_L[bit_L, :]), @view(f_vs_data.midpoint[offset+1:end, :]))
    f_vs_data.vn = @view(f_vs_data.midpoint[:, DIR])
    f_vs_data.df = Vcat(df_L, @view(f_vs_data.df[offset+1:end, :]))
    f_vs_data.sdf =
        Vcat(@view(vs_data.sdf[bit_L, :, DIR]), @view(f_vs_data.sdf[offset+1:end, :, DIR]))
    f_vs_data.weight =
        Vcat(@view(vs_data.weight[bit_L]), @view(f_vs_data.weight[offset+1:end]))
    offset = length(df_L)
    return (bit_L, offset)
end

function calc_flux!(::SameSizeNeighbor, face::Face{InnerFace}, amr::AMR{2,2})
    ps_data = face.data
    faceid = face.faceid
    DIR = get_dir(faceid)
    ROT = get_rot(faceid) # Left: -1, Right: 1
    global_data = amr.global_data
    gas = global_data.config.gas
    ds = ps_data.ds[DIR]
    vs_data = ps_data.vs_data
    nps_data = ps_data.neighbor.data[faceid][1]
    vs_data_n = nps_data.vs_data
    f_vs_data, offset, bit_L, _ = make_face_data(ps_data, nps_data, faceid, ROT, DIR)
    reconstruct_vs!(f_vs_data, ds, ds, offset, ROT)
    w0 = calc_w0(f_vs_data)
    prim0 = get_prim(w0, global_data)
    qf0 = calc_qf(f_vs_data, prim0)
    aL, aR = calc_a(w0, prim0, ps_data.w, nps_data.w, ds, ds, global_data, ROT)
    Mu, Mv, Mξ, Mu_L, Mu_R = moment_u(prim0, global_data, ROT, DIR)
    A = calc_A(prim0, aL, aR, Mu, Mv, Mξ, Mu_L, Mu_R, global_data, DIR)
    τ0 = get_τ(prim0, gas.μᵣ, gas.ω)
    Mt = time_int(τ0, global_data.status.Δt)
    fw = calc_flux_g0_2D2F(prim0, Mt, Mu, Mv, Mξ, Mu_L, Mu_R, aL, aR, A, DIR)
    F = discrete_maxwell(f_vs_data.midpoint, prim0, global_data)
    F⁺ = shakhov_part(f_vs_data.midpoint, F, prim0, qf0, global_data)
    fw .+= calc_flux_f0(f_vs_data, F⁺, Mt)
    dsf = calc_face_area(ps_data, DIR)
    fw .*= dsf * ROT
    ps_data.flux .+= fw
    update_nflux!(nps_data, fw)
    micro_flux = calc_micro_flux(f_vs_data, F, F⁺, aL, aR, A, Mξ, Mt, offset, dsf)
    update_vs_flux!(micro_flux, bit_L, vs_data, vs_data_n, offset, ROT)
end
function calc_flux!(::SameSizeNeighbor, face::Face{InnerFace}, amr::AMR{3,1})
    ps_data = face.data
    faceid = face.faceid
    DIR = get_dir(faceid)
    ROT = get_rot(faceid) # Left: -1, Right: 1
    global_data = amr.global_data
    gas = global_data.config.gas
    ds = ps_data.ds[DIR]
    vs_data = ps_data.vs_data
    nps_data = ps_data.neighbor.data[faceid][1]
    vs_data_n = nps_data.vs_data
    f_vs_data, offset, bit_L, _ = make_face_data(ps_data, nps_data, faceid, ROT, DIR)
    reconstruct_vs!(f_vs_data, ds, ds, offset, ROT)
    w0 = calc_w0(f_vs_data)
    prim0 = get_prim(w0, global_data)
    qf0 = calc_qf(f_vs_data, prim0)
    aL, aR = calc_a(w0, prim0, ps_data.w, nps_data.w, ds, ds, global_data, ROT)
    Mu, Mv, Mw, Mu_L, Mu_R = moment_u(prim0, global_data, ROT, DIR)
    A = calc_A(prim0, aL, aR, Mu, Mv, Mw, Mu_L, Mu_R, global_data, DIR)
    τ0 = get_τ(prim0, gas.μᵣ, gas.ω)
    Mt = time_int(τ0, global_data.status.Δt)
    fw = calc_flux_g0_3D1F(prim0, Mt, Mu, Mv, Mw, Mu_L, Mu_R, aL, aR, A, DIR)
    F = discrete_maxwell(f_vs_data.midpoint, prim0, global_data)
    F⁺ = shakhov_part(f_vs_data.midpoint, F, prim0, qf0, global_data)
    fw .+= calc_flux_f0(f_vs_data, F⁺, Mt)
    dsf = calc_face_area(ps_data, DIR)
    fw .*= dsf * ROT
    ps_data.flux .+= fw
    update_nflux!(nps_data, fw)
    micro_flux = calc_micro_flux(f_vs_data, F, F⁺, aL, aR, A, Mt, offset, dsf)
    update_vs_flux!(micro_flux, bit_L, vs_data, vs_data_n, offset, ROT)
end
function calc_flux!(::HalfSizeNeighbor, face::Face, amr::AMR{2,2})
    ps_data = face.data
    faceid = face.faceid
    DIR = get_dir(faceid)
    ROT = get_rot(faceid) # Left: -1, Right: 1
    global_data = amr.global_data
    gas = global_data.config.gas
    ds = ps_data.ds[DIR]
    vs_data = ps_data.vs_data
    for i = 1:2
        nps_data = ps_data.neighbor.data[faceid][i]
        vs_data_n = nps_data.vs_data
        if i == 1
            f_vs_data, offset, bit_L, _ =
                make_face_data(ps_data, nps_data, faceid, ROT, DIR)
            reconstruct_vs!(f_vs_data, ds, 0.5 * ds, offset, ROT)
        else
            update_face_data_R!(ps_data, nps_data, f_vs_data, offset, faceid, ROT, DIR)
            reconstruct_vs_R!(f_vs_data, 0.5 * ds, offset, ROT)
        end
        w0 = calc_w0(f_vs_data)
        prim0 = get_prim(w0, global_data)
        qf0 = calc_qf(f_vs_data, prim0)
        aL, aR = calc_a(w0, prim0, ps_data.w, nps_data.w, ds, 0.5 * ds, global_data, ROT)
        Mu, Mv, Mξ, Mu_L, Mu_R = moment_u(prim0, global_data, ROT, DIR)
        A = calc_A(prim0, aL, aR, Mu, Mv, Mξ, Mu_L, Mu_R, global_data, DIR)
        τ0 = get_τ(prim0, gas.μᵣ, gas.ω)
        Mt = time_int(τ0, global_data.status.Δt)
        fw = calc_flux_g0_2D2F(prim0, Mt, Mu, Mv, Mξ, Mu_L, Mu_R, aL, aR, A, DIR)
        F = discrete_maxwell(f_vs_data.midpoint, prim0, global_data)
        F⁺ = shakhov_part(f_vs_data.midpoint, F, prim0, qf0, global_data)
        fw .+= calc_flux_f0(f_vs_data, F⁺, Mt)
        dsf = calc_face_area(nps_data, DIR)
        fw .*= dsf * ROT
        ps_data.flux .+= fw
        update_nflux!(nps_data, fw)
        micro_flux = calc_micro_flux(f_vs_data, F, F⁺, aL, aR, A, Mξ, Mt, offset, dsf)
        update_vs_flux!(micro_flux, bit_L, vs_data, vs_data_n, offset, ROT)
    end
end
function calc_flux!(::HalfSizeNeighbor, face::Face, amr::AMR{3,1})
    ps_data = face.data
    faceid = face.faceid
    DIR = get_dir(faceid)
    ROT = get_rot(faceid) # Left: -1, Right: 1
    global_data = amr.global_data
    gas = global_data.config.gas
    ds = ps_data.ds[DIR]
    vs_data = ps_data.vs_data
    f_vs_data = nothing; offset = 0; bit_L = nothing
    for i = 1:4
        nps_data = ps_data.neighbor.data[faceid][i]
        vs_data_n = nps_data.vs_data
        if i == 1
            f_vs_data, offset, bit_L, _ =
                make_face_data(ps_data, nps_data, faceid, ROT, DIR)
            reconstruct_vs!(f_vs_data, ds, 0.5 * ds, offset, ROT)
        else
            update_face_data_R!(ps_data, nps_data, f_vs_data, offset, faceid, ROT, DIR)
            reconstruct_vs_R!(f_vs_data, 0.5 * ds, offset, ROT)
        end
        w0 = calc_w0(f_vs_data)
        prim0 = get_prim(w0, global_data)
        qf0 = calc_qf(f_vs_data, prim0)
        aL, aR = calc_a(w0, prim0, ps_data.w, nps_data.w, ds, 0.5 * ds, global_data, ROT)
        Mu, Mv, Mw, Mu_L, Mu_R = moment_u(prim0, global_data, ROT, DIR)
        A = calc_A(prim0, aL, aR, Mu, Mv, Mw, Mu_L, Mu_R, global_data, DIR)
        τ0 = get_τ(prim0, gas.μᵣ, gas.ω)
        Mt = time_int(τ0, global_data.status.Δt)
        fw = calc_flux_g0_3D1F(prim0, Mt, Mu, Mv, Mw, Mu_L, Mu_R, aL, aR, A, DIR)
        F = discrete_maxwell(f_vs_data.midpoint, prim0, global_data)
        F⁺ = shakhov_part(f_vs_data.midpoint, F, prim0, qf0, global_data)
        fw .+= calc_flux_f0(f_vs_data, F⁺, Mt)
        dsf = calc_face_area(nps_data, DIR)
        fw .*= dsf * ROT
        ps_data.flux .+= fw
        update_nflux!(nps_data, fw)
        micro_flux = calc_micro_flux(f_vs_data, F, F⁺, aL, aR, A, Mt, offset, dsf)
        update_vs_flux!(micro_flux, bit_L, vs_data, vs_data_n, offset, ROT)
    end
end
function calc_flux!(::DoubleSizeNeighbor, face::Face, amr::AMR{2,2})
    ps_data = face.data
    faceid = face.faceid
    DIR = get_dir(faceid)
    ROT = get_rot(faceid) # Left: -1, Right: 1
    global_data = amr.global_data
    gas = global_data.config.gas
    ds = ps_data.ds[DIR]
    vs_data = ps_data.vs_data
    nps_data = ps_data.neighbor.data[faceid][1]
    vs_data_n = nps_data.vs_data
    f_vs_data, offset, bit_L, _ = make_face_data(ps_data, nps_data, faceid, ROT, DIR)
    reconstruct_vs!(f_vs_data, ds, 2.0 * ds, offset, ROT)
    w0 = calc_w0(f_vs_data)
    prim0 = get_prim(w0, global_data)
    qf0 = calc_qf(f_vs_data, prim0)
    aL, aR = calc_a(w0, prim0, ps_data.w, nps_data.w, ds, 2.0 * ds, global_data, ROT)
    Mu, Mv, Mw, Mu_L, Mu_R = moment_u(prim0, global_data, ROT, DIR)
    A = calc_A(prim0, aL, aR, Mu, Mv, Mw, Mu_L, Mu_R, global_data, DIR)
    τ0 = get_τ(prim0, gas.μᵣ, gas.ω)
    Mt = time_int(τ0, global_data.status.Δt)
    fw = calc_flux_g0_2D2F(prim0, Mt, Mu, Mv, Mw, Mu_L, Mu_R, aL, aR, A, DIR)
    F = discrete_maxwell(f_vs_data.midpoint, prim0, global_data)
    F⁺ = shakhov_part(f_vs_data.midpoint, F, prim0, qf0, global_data)
    fw .+= calc_flux_f0(f_vs_data, F⁺, Mt)
    dsf = calc_face_area(ps_data, DIR)
    fw .*= dsf * ROT
    ps_data.flux .+= fw
    micro_flux = calc_micro_flux(f_vs_data, F, F⁺, aL, aR, A, Mt, offset, dsf)
    update_vs_flux!(micro_flux, bit_L, vs_data, vs_data_n, offset, ROT)
    for i in eachindex(face.hanging_data)
        (isa(face.hanging_data[i], Ghost_PS_Data)||isa(face.hanging_data[i],MissingHangingQuad)) && continue
        ps_data = face.hanging_data[i]
        bit_L, offset = update_face_data_L!(ps_data, nps_data, f_vs_data, offset, faceid, ROT, DIR)
        reconstruct_vs_L!(f_vs_data, ds, offset, ROT)
        w0 = calc_w0(f_vs_data)
        prim0 = get_prim(w0, global_data)
        qf0 = calc_qf(f_vs_data, prim0)
        aL, aR = calc_a(w0, prim0, ps_data.w, nps_data.w, ds, 2.0 * ds, global_data, ROT)
        Mu, Mv, Mξ, Mu_L, Mu_R = moment_u(prim0, global_data, ROT, DIR)
        A = calc_A(prim0, aL, aR, Mu, Mv, Mξ, Mu_L, Mu_R, global_data, DIR)
        τ0 = get_τ(prim0, gas.μᵣ, gas.ω)
        Mt = time_int(τ0, global_data.status.Δt)
        fw = calc_flux_g0_2D2F(prim0, Mt, Mu, Mv, Mξ, Mu_L, Mu_R, aL, aR, A, DIR)
        F = discrete_maxwell(f_vs_data.midpoint, prim0, global_data)
        F⁺ = shakhov_part(f_vs_data.midpoint, F, prim0, qf0, global_data)
        fw .+= calc_flux_f0(f_vs_data, F⁺, Mt)
        fw .*= dsf * ROT
        ps_data.flux .+= fw
        micro_flux = calc_micro_flux(f_vs_data, F, F⁺, aL, aR, A, Mξ, Mt, offset, dsf)
    end
end
function calc_flux!(::DoubleSizeNeighbor, face::Face, amr::AMR{3,1})
    ps_data = face.data
    faceid = face.faceid
    DIR = get_dir(faceid)
    ROT = get_rot(faceid) # Left: -1, Right: 1
    global_data = amr.global_data
    gas = global_data.config.gas
    ds = ps_data.ds[DIR]
    vs_data = ps_data.vs_data
    nps_data = ps_data.neighbor.data[faceid][1]
    vs_data_n = nps_data.vs_data
    f_vs_data, offset, bit_L, _ = make_face_data(ps_data, nps_data, faceid, ROT, DIR)
    reconstruct_vs!(f_vs_data, ds, 2.0 * ds, offset, ROT)
    w0 = calc_w0(f_vs_data)
    prim0 = get_prim(w0, global_data)
    qf0 = calc_qf(f_vs_data, prim0)
    aL, aR = calc_a(w0, prim0, ps_data.w, nps_data.w, ds, 2.0 * ds, global_data, ROT)
    Mu, Mv, Mw, Mu_L, Mu_R = moment_u(prim0, global_data, ROT, DIR)
    A = calc_A(prim0, aL, aR, Mu, Mv, Mw, Mu_L, Mu_R, global_data, DIR)
    τ0 = get_τ(prim0, gas.μᵣ, gas.ω)
    Mt = time_int(τ0, global_data.status.Δt)
    fw = calc_flux_g0_3D1F(prim0, Mt, Mu, Mv, Mw, Mu_L, Mu_R, aL, aR, A, DIR)
    F = discrete_maxwell(f_vs_data.midpoint, prim0, global_data)
    F⁺ = shakhov_part(f_vs_data.midpoint, F, prim0, qf0, global_data)
    fw .+= calc_flux_f0(f_vs_data, F⁺, Mt)
    dsf = calc_face_area(ps_data, DIR)
    fw .*= dsf * ROT
    ps_data.flux .+= fw
    micro_flux = calc_micro_flux(f_vs_data, F, F⁺, aL, aR, A, Mt, offset, dsf)
    update_vs_flux!(micro_flux, bit_L, vs_data, vs_data_n, offset, ROT)
    for i in eachindex(face.hanging_data)
        (isa(face.hanging_data[i], Ghost_PS_Data)||isa(face.hanging_data[i],MissingHangingQuad)) && continue
        ps_data = face.hanging_data[i]
        bit_L, offset = update_face_data_L!(ps_data, nps_data, f_vs_data, offset, faceid, ROT, DIR)
        reconstruct_vs_L!(f_vs_data, ds, offset, ROT)
        w0 = calc_w0(f_vs_data)
        prim0 = get_prim(w0, global_data)
        qf0 = calc_qf(f_vs_data, prim0)
        aL, aR = calc_a(w0, prim0, ps_data.w, nps_data.w, ds, 2.0 * ds, global_data, ROT)
        Mu, Mv, Mξ, Mu_L, Mu_R = moment_u(prim0, global_data, ROT, DIR)
        A = calc_A(prim0, aL, aR, Mu, Mv, Mw, Mu_L, Mu_R, global_data, DIR)
        τ0 = get_τ(prim0, gas.μᵣ, gas.ω)
        Mt = time_int(τ0, global_data.status.Δt)
        fw = calc_flux_g0_3D1F(prim0, Mt, Mu, Mv, Mw, Mu_L, Mu_R, aL, aR, A, DIR)
        F = discrete_maxwell(f_vs_data.midpoint, prim0, global_data)
        F⁺ = shakhov_part(f_vs_data.midpoint, F, prim0, qf0, global_data)
        fw .+= calc_flux_f0(f_vs_data, F⁺, Mt)
        fw .*= dsf * ROT
        ps_data.flux .+= fw
        micro_flux = calc_micro_flux(f_vs_data, F, F⁺, aL, aR, A, Mt, offset, dsf)
    end
end
function update_flux!(amr::AMR)
    faces = amr.field.faces
    @inbounds @simd for i in eachindex(faces)
        face = faces[i]
        calc_flux!(Val(face.data.neighbor.state[face.faceid]), face, amr)
    end
end