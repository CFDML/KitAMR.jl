function calc_flux!(::Val{1}, ::Val{0}, face::Face, DVM_data::DVM_Data, faceid::Int)
    ps_data = face.data
    dir = get_dir(faceid)
    rot = get_rot(faceid)
    vs_data = ps_data.vs_data
    midpoint = vs_data.midpoint
    gas = DVM_data.global_data.gas
    Θ = [heaviside(rot * midpoint[i, dir]) for i in axes(midpoint, 1)]
    h = @views vs_data.df[:, 1]
    u = @views midpoint[:, 1]
    v = @views midpoint[:, 2]
    vn = @views midpoint[:, dir]
    prim0 = DVM_data.global_data.bc[:, faceid]
    df0 = @. vs_data.df - 0.5 * rot * ps_data.ds[dir] * @view(vs_data.sdf[:, :, dir])
    h0 = @views df0[:, 1]
    b0 = @views df0[:, 2]
    prim0[1] = calc_ρw(prim0, h0, Θ, u, v, vs_data.weight, vn)
    H0, B0 = discrete_maxwell(u, v, prim0, gas.K)
    H = @. H0 * Θ + h0 * (1.0 - Θ)
    B = @. B0 * Θ + b0 * (1.0 - Θ)
    fw = calc_fwb(H, B, u, v, vs_data.weight, vn)
    @. ps_data.flux += rot * fw * ps_data.ds[dir] * gas.Δt
    @. vs_data.flux[:, 1] += vn * rot * ps_data.ds[dir] * gas.Δt * H
    # @. vs_data.flux[:,1] += 1.
    # if ps_data.midpoint==[0.011111111111111112, 0.9888888888888889]
    #     @show findmax(vs_data.flux[:,1]) faceid
    # end
    @. vs_data.flux[:, 2] += vn * rot * ps_data.ds[dir] * gas.Δt * B
end
function make_face_data(ps_data::PS_Data, nps_data::AbstractPsData, rot::Float64, dir::Int)
    vs_data = ps_data.vs_data
    midpoint_L = vs_data.midpoint
    vs_data_n = nps_data.vs_data
    midpoint_R = vs_data_n.midpoint
    bit_L = [x <= 0.0 for x in rot .* @view(midpoint_L[:, dir])]
    bit_R = [x > 0.0 for x in rot .* @view(midpoint_R[:, dir])]
    midpoint_n = vcat(midpoint_L[bit_L, :], midpoint_R[bit_R, :])
    df_n = vcat(vs_data.df[bit_L, :], vs_data_n.df[bit_R, :])
    sdf_n = vcat(vs_data.sdf[bit_L, :, dir], vs_data_n.sdf[bit_R, :, dir])
    weight_L = vs_data.weight[bit_L]
    offset = length(weight_L)
    weight_n = vcat(weight_L, vs_data_n.weight[bit_R])
    (weight_n, midpoint_n, df_n, sdf_n, offset, bit_L, bit_R)
end
# function update_nvs_flux!(::Val{1},fh::AV,fb,::Int,vs_data_n,rot::Float64,nid::Vector{Int},index::Int)
#     vs_data_n.flux[nid[2],1] -= rot*fh[index]
#     vs_data_n.flux[nid[2],2] -= rot*fb[index]
#     nid[2] += 1
# end
# function update_nvs_flux!(::Val{2^DIM},fh::AV,fb,::Int,vs_data_n,rot::Float64,nid::Vector{Int},index::Int)
#     for _ = 1:2^DIM
#         vs_data_n.flux[nid[2],1] -= rot*fh[index]/2^DIM
#         vs_data_n.flux[nid[2],2] -= rot*fb[index]/2^DIM
#         nid[2] += 1
#     end
# end
# function update_nvs_flux!(::Val{-1},fh::AV,fb,lindex::Int,vs_data_n,rot::Float64,nid::Vector{Int},index::Int)
#     vs_data_n.flux[nid[2],1] -= rot*fh[index]
#     vs_data_n.flux[nid[2],2] -= rot*fb[index]
#    (lindex===2^DIM) && (nid[2] += 1)
# end
# function update_vs_flux!(fh::AV,fb::AV,bit_L::AV,bit_R::AV,vs_data::VS_Data,vs_data_n::VS_Data,state,state_n,lindex::Vector{Int},nindex::Vector{Int},rot::Float64,)
#     nid = [1,1]
#     @inbounds @simd for i in eachindex(bit_L)
#         if bit_L[i]
#             index = nid[1]
#             vs_data.flux[i,1] += rot*fh[index]
#             vs_data.flux[i,2] += rot*fb[index]
#             update_nvs_flux!(Val(state[i]),fh,fb,lindex[i],vs_data_n,rot,nid,index)
#             nid[1] += 1
#         else
#             nid[2] += state[i]
#         end
#     end
#     nid[2] = 1
#     rot = -rot
#     @inbounds @simd for i in eachindex(bit_R)
#         if bit_R[i]
#             index = nid[1]
#             vs_data_n.flux[i,1] += rot*fh[index]
#             vs_data_n.flux[i,2] += rot*fb[index]
#             update_nvs_flux!(Val(state_n[i]),fh,fb,nindex[i],vs_data,rot,nid,index)
#             nid[1] += 1
#         else
#             nid[2] += state_n[i]
#         end
#     end
# end
# function update_vs_flux!(fh::AV,fb::AV,bit_L::AV,bit_R::AV,vs_data::VS_Data,::Ghost_VS_Data,state,state_n,lindex::Vector{Int},nindex::Vector{Int},rot::Float64,)
#     nid = [1,1]
#     @inbounds @simd for i in eachindex(bit_L)
#         if bit_L[i]
#             index = nid[1]
#             vs_data.flux[i,1] += rot*fh[index]
#             vs_data.flux[i,2] += rot*fb[index]
#             nid[1] += 1
#         else
#             nid[2] += state[i]
#         end
#     end
#     nid[2] = 1
#     rot = -rot
#     @inbounds @simd for i in eachindex(bit_R)
#         if bit_R[i]
#             index = nid[1]
#             # state = vs_structure_n[i]
#             if state_n[i]>1
#                 @show state_n
#             end
#             update_nvs_flux!(Val(state_n[i]),fh,fb,nindex[i],vs_data,rot,nid,index)
#             nid[1] += 1
#         else
#             nid[2] += state_n[i]
#         end
#     end
# end
function update_vs_flux!(
    fh::AV,
    fb::AV,
    bit_L::AV,
    vs_data::AbstractVsData,
    vs_data_n::VS_Data,
    offset::Int,
    rot::Float64,
)
    index = j = 1
    index_n = offset + 1
    flag = 0.0
    level = vs_data.level
    level_n = vs_data_n.level
    flux = vs_data.flux
    flux_n = vs_data_n.flux
    for i = 1:vs_data.vs_num
        if bit_L[i]
            flux[i, 1] += rot * fh[index]
            flux[i, 2] += rot * fb[index]
            if level[i] == level_n[j]
                flux_n[j, 1] -= rot * fh[index]
                flux_n[j, 2] -= rot * fb[index]
                j += 1
            elseif level[i] < level_n[j]
                while flag != 1.0
                    flux_n[j, 1] -= rot * fh[index]
                    flux_n[j, 2] -= rot * fb[index]
                    flag += 1 / 2^(DIM + level_n[j] - level[i])
                    j += 1
                end
                flag = 0.0
            else
                flux_n[j, 1] -= rot * fh[index] / 2^DIM
                flux_n[j, 2] -= rot * fb[index] / 2^DIM
                flag += 1 / 2^(DIM + level[i] - level_n[j])
                if flag == 1.0
                    j += 1
                    flag = 0.0
                end
            end
            index += 1
        else
            if level[i] == level_n[j]
                flux_n[j, 1] -= rot * fh[index_n]
                flux_n[j, 2] -= rot * fb[index_n]
                flux[i, 1] += rot * fh[index_n]
                flux[i, 2] += rot * fb[index_n]
                j += 1
                index_n += 1
            elseif level[i] < level_n[j]
                while flag != 1.0
                    flux_n[j, 1] -= rot * fh[index_n]
                    flux_n[j, 2] -= rot * fb[index_n]
                    flux[i, 1] += rot * fh[index_n] / 2^DIM
                    flux[i, 2] += rot * fb[index_n] / 2^DIM
                    flag += 1 / 2^(DIM + level_n[j] - level[i])
                    j += 1
                    index_n += 1
                end
                flag = 0.0
            else
                flux_n[j, 1] -= rot * fh[index_n] / 2^DIM
                flux_n[j, 2] -= rot * fb[index_n] / 2^DIM
                flux[i, 1] += rot * fh[index_n]
                flux[i, 2] += rot * fb[index_n]
                flag += 1 / 2^(DIM + level[i] - level_n[index])
                if flag == 1.0
                    j += 1
                    index_n += 1
                    flag = 0.0
                end
            end
        end
    end
end
function update_vs_flux!(
    fh::AV,
    fb::AV,
    bit_L::AV,
    vs_data::AbstractVsData,
    vs_data_n::Ghost_VS_Data,
    offset::Int,
    rot::Float64,
)
    index = j = 1
    index_n = offset + 1
    flag = 0.0
    level = vs_data.level
    level_n = vs_data_n.level
    flux = vs_data.flux
    for i = 1:vs_data.vs_num
        if bit_L[i]
            flux[i, 1] += rot * fh[index]
            flux[i, 2] += rot * fb[index]
            if level[i] == level_n[j]
                j += 1
            elseif level[i] < level_n[j]
                while flag != 1.0
                    flag += 1 / 2^(DIM + level_n[j] - level[i])
                    j += 1
                end
                flag = 0.0
            else
                flag += 1 / 2^(DIM + level[i] - level_n[j])
                if flag == 1.0
                    j += 1
                    flag = 0.0
                end
            end
            index += 1
        else
            if level[i] == level_n[j]
                flux[i, 1] += rot * fh[index_n]
                flux[i, 2] += rot * fb[index_n]
                j += 1
                index_n += 1
            elseif level[i] < level_n[j]
                while flag != 1.0
                    flux[i, 1] += rot * fh[index_n] / 2^DIM
                    flux[i, 2] += rot * fb[index_n] / 2^DIM
                    flag += 1 / 2^(DIM + level_n[j] - level[i])
                    j += 1
                    index_n += 1
                end
                flag = 0.0
            else
                flux[i, 1] += rot * fh[index_n]
                flux[i, 2] += rot * fb[index_n]
                flag += 1 / 2^(DIM + level[i] - level_n[index])
                if flag == 1.0
                    j += 1
                    index_n += 1
                    flag = 0.0
                end
            end
        end
    end
end
function update_nflux!(nps_data::PS_Data, fw::AV)
    nps_data.flux .-= fw
end
function update_nflux!(::Ghost_PS_Data, ::AV)
    return nothing
end
function calc_flux!(::Val{0}, ::Val{1}, face::Face, DVM_data::DVM_Data, faceid::Int)
    ps_data = face.data
    dir = get_dir(faceid)
    rot = get_rot(faceid) # Left: -1, Right: 1
    gas = DVM_data.global_data.gas
    ds = ps_data.ds[dir]
    vs_data = ps_data.vs_data
    nps_data = ps_data.neighbor.data[faceid][1]
    vs_data_n = nps_data.vs_data
    weight_n, midpoint_n, df_n, sdf_n, offset, bit_L, _ =
        make_face_data(ps_data, nps_data, rot, dir)
    h = @view(df_n[:, 1])
    b = @view(df_n[:, 2])
    u = @view(midpoint_n[:, 1])
    v = @view(midpoint_n[:, 2])
    vn = @view(midpoint_n[:, dir])
    sh = @view(sdf_n[:, 1])
    sb = @view(sdf_n[:, 2])
    h0, b0 = reconstruct_vs(h, b, sh, sb, ds, ds, offset, rot)
    w0 = calc_w0(h0, b0, u, v, weight_n)
    prim0 = get_prim(w0, gas.γ)
    # if ps_data.midpoint==[0.011111111111111112, 0.9888888888888889]&&faceid==2
    #     @show w0 prim0 maximum(vs_data.sdf[:,1,1]) maximum(vs_data.df[:,1]) maximum(vs_data_n.df[:,1])
    # end
    qf0 = calc_qf(h0, b0, prim0, u, v, weight_n)
    aL, aR = calc_a(w0, prim0, ps_data.w, nps_data.w, ds, ds, gas.K, rot)
    Mu, Mv, Mξ, Mu_L, Mu_R = moment_u(prim0, 6, 4, gas.K, rot, dir)
    A = calc_A(prim0, aL, aR, Mu, Mv, Mu_L, Mu_R, Mξ, gas.K, dir)
    τ0 = get_τ(prim0, gas.μᵣ, gas.ω)
    Mt = calc_time_int(τ0, gas.Δt)
    fw = calc_flux_g0(prim0, Mt, Mu, Mu_L, Mu_R, Mv, Mξ, aL, aR, A, dir)
    H, B = discrete_maxwell(u, v, prim0, gas.K)
    H⁺, B⁺ = shakhov_part(H, B, prim0, u, v, qf0, gas.Pr, gas.K)
    fw .+= calc_flux_f0(Mt, H⁺, B⁺, u, v, weight_n, h0, b0, sh, sb, vn)
    fw .*= ds * rot
    ps_data.flux .+= fw
    update_nflux!(nps_data, fw)
    fh, fb = calc_fhb(H, B, H⁺, B⁺, Mt, Mξ, u, v, h0, b0, sh, sb, aL, aR, A, offset, vn, ds)
    # if ps_data.midpoint==[0.011111111111111112, 0.9888888888888889]&&faceid==2
    #     @show maximum(fh) maximum(vs_data.flux[:,1])
    # #    @show maximum(vs_data.flux[:,1])
    # end
    # fh_old = vs_data.flux[:,1]
    # fh = ones(vs_data.vs_num)
    # fh_old = vs_data_n.flux[:,1]
    # if nps_data.midpoint == [0.90625, 0.90625]
    #     @show nps_data.flux faceid
    # end
    update_vs_flux!(fh, fb, bit_L, vs_data, vs_data_n, offset, rot)
    # if nps_data.midpoint==[0.011111111111111112, 0.9888888888888889]
    #     # @show aL aR Mu_L Mt A f1 f2 fw
    #     @show findmax(vs_data_n.flux[:,1]) faceid findmax(fh) bit_L[436] fh[436] vs_data_n.flux[436,1]+fh[436]
    # end
    # if ps_data.midpoint==[0.011111111111111112, 0.9888888888888889]
    #     @show findmax(vs_data.flux[:,1]) faceid
    #     # @show maximum(vs_data.sdf[:,1,1]) maximum(vs_data.sdf[:,1,2])
    # end
end
function calc_flux!(::Val{0}, ::Val{2}, face::Face, DVM_data::DVM_Data, faceid::Int)
    ps_data = face.data
    dir = get_dir(faceid)
    rot = get_rot(faceid) # Left: -1, Right: 1
    gas = DVM_data.global_data.gas
    vs_data = ps_data.vs_data
    # if ps_data.midpoint == [0.90625, 0.90625]
    #     @show ps_data.neighbor.data[faceid][1].midpoint
    #     @show ps_data.neighbor.data[faceid][2].midpoint
    # end
    fw_test = zeros(DIM + 2)
    for i = 1:2^(DIM-1)
        nps_data = ps_data.neighbor.data[faceid][i]
        vs_data_n = nps_data.vs_data
        ds = nps_data.ds[dir]
        weight_n, midpoint_n, df_n, sdf_n, offset, bit_L, _ =
            make_face_data(ps_data, nps_data, rot, dir)
        h = @view(df_n[:, 1])
        b = @view(df_n[:, 2])
        u = @view(midpoint_n[:, 1])
        v = @view(midpoint_n[:, 2])
        vn = @view(midpoint_n[:, dir])
        sh = @view(sdf_n[:, 1])
        sb = @view(sdf_n[:, 2])
        h0, b0 = reconstruct_vs(h, b, sh, sb, 2.0 * ds, ds, offset, rot)
        w0 = calc_w0(h0, b0, u, v, weight_n)
        prim0 = get_prim(w0, gas.γ)
        qf0 = calc_qf(h0, b0, prim0, u, v, weight_n)
        aL, aR = calc_a(w0, prim0, ps_data.w, nps_data.w, 2.0 * ds, ds, gas.K, rot)
        Mu, Mv, Mξ, Mu_L, Mu_R = moment_u(prim0, 6, 4, gas.K, rot, dir)
        A = calc_A(prim0, aL, aR, Mu, Mv, Mu_L, Mu_R, Mξ, gas.K, dir)
        τ0 = get_τ(prim0, gas.μᵣ, gas.ω)
        Mt = calc_time_int(τ0, gas.Δt)
        fw = calc_flux_g0(prim0, Mt, Mu, Mu_L, Mu_R, Mv, Mξ, aL, aR, A, dir)
        H, B = discrete_maxwell(u, v, prim0, gas.K)
        H⁺, B⁺ = shakhov_part(H, B, prim0, u, v, qf0, gas.Pr, gas.K)
        fw .+= calc_flux_f0(Mt, H⁺, B⁺, u, v, weight_n, h0, b0, sh, sb, vn)
        fw .*= ds * rot
        ps_data.flux .+= fw
        fw_test .+= fw
        update_nflux!(nps_data, fw)
        fh, fb =
            calc_fhb(H, B, H⁺, B⁺, Mt, Mξ, u, v, h0, b0, sh, sb, aL, aR, A, offset, vn, ds)
        update_vs_flux!(fh, fb, bit_L, vs_data, vs_data_n, offset, rot)
    end
    # if ps_data.midpoint == [0.90625, 0.90625]
    #     @show ps_data.flux faceid
    # end
end
function calc_flux!(::Val{2}, ::Val{-1}, face::Face, DVM_data::DVM_Data, faceid::Int)
    ps_data = face.data
    dir = get_dir(faceid)
    rot = get_rot(faceid) # Left: -1, Right: 1
    gas = DVM_data.global_data.gas
    ds = ps_data.ds[dir]
    vs_data = ps_data.vs_data
    nps_data = ps_data.neighbor.data[faceid][1]
    vs_data_n = nps_data.vs_data
    weight_n, midpoint_n, df_n, sdf_n, offset, bit_L, _ =
        make_face_data(ps_data, nps_data, rot, dir)
    h = @view(df_n[:, 1])
    b = @view(df_n[:, 2])
    u = @view(midpoint_n[:, 1])
    v = @view(midpoint_n[:, 2])
    vn = @view(midpoint_n[:, dir])
    sh = @view(sdf_n[:, 1])
    sb = @view(sdf_n[:, 2])
    h0, b0 = reconstruct_vs(h, b, sh, sb, ds, 2.0 * ds, offset, rot)
    w0 = calc_w0(h0, b0, u, v, weight_n)
    prim0 = get_prim(w0, gas.γ)
    qf0 = calc_qf(h0, b0, prim0, u, v, weight_n)
    aL, aR = calc_a(w0, prim0, ps_data.w, nps_data.w, ds, 2.0 * ds, gas.K, rot)
    Mu, Mv, Mξ, Mu_L, Mu_R = moment_u(prim0, 6, 4, gas.K, rot, dir)
    A = calc_A(prim0, aL, aR, Mu, Mv, Mu_L, Mu_R, Mξ, gas.K, dir)
    τ0 = get_τ(prim0, gas.μᵣ, gas.ω)
    Mt = calc_time_int(τ0, gas.Δt)
    fw = calc_flux_g0(prim0, Mt, Mu, Mu_L, Mu_R, Mv, Mξ, aL, aR, A, dir)
    H, B = discrete_maxwell(u, v, prim0, gas.K)
    H⁺, B⁺ = shakhov_part(H, B, prim0, u, v, qf0, gas.Pr, gas.K)
    fw .+= calc_flux_f0(Mt, H⁺, B⁺, u, v, weight_n, h0, b0, sh, sb, vn)
    fw .*= ds * rot
    ps_data.flux .+= fw
    fh, fb = calc_fhb(H, B, H⁺, B⁺, Mt, Mξ, u, v, h0, b0, sh, sb, aL, aR, A, offset, vn, ds)
    update_vs_flux!(fh, fb, bit_L, vs_data, vs_data_n, offset, rot)
    ps_data = face.hanging_data
    vs_data = ps_data.vs_data
    weight_n, midpoint_n, df_n, sdf_n, offset, bit_L, _ =
        make_face_data(ps_data, nps_data, rot, dir)
    h = @view(df_n[:, 1])
    b = @view(df_n[:, 2])
    u = @view(midpoint_n[:, 1])
    v = @view(midpoint_n[:, 2])
    vn = @view(midpoint_n[:, dir])
    sh = @view(sdf_n[:, 1])
    sb = @view(sdf_n[:, 2])
    h0, b0 = reconstruct_vs(h, b, sh, sb, ds, 2.0 * ds, offset, rot)
    w0 = calc_w0(h0, b0, u, v, weight_n)
    prim0 = get_prim(w0, gas.γ)
    qf0 = calc_qf(h0, b0, prim0, u, v, weight_n)
    aL, aR = calc_a(w0, prim0, ps_data.w, nps_data.w, ds, 2.0 * ds, gas.K, rot)
    Mu, Mv, Mξ, Mu_L, Mu_R = moment_u(prim0, 6, 4, gas.K, rot, dir)
    A = calc_A(prim0, aL, aR, Mu, Mv, Mu_L, Mu_R, Mξ, gas.K, dir)
    τ0 = get_τ(prim0, gas.μᵣ, gas.ω)
    Mt = calc_time_int(τ0, gas.Δt)
    fw = calc_flux_g0(prim0, Mt, Mu, Mu_L, Mu_R, Mv, Mξ, aL, aR, A, dir)
    H, B = discrete_maxwell(u, v, prim0, gas.K)
    H⁺, B⁺ = shakhov_part(H, B, prim0, u, v, qf0, gas.Pr, gas.K)
    fw .+= calc_flux_f0(Mt, H⁺, B⁺, u, v, weight_n, h0, b0, sh, sb, vn)
    fw .*= ds * rot
    ps_data.flux .+= fw
    fh, fb = calc_fhb(H, B, H⁺, B⁺, Mt, Mξ, u, v, h0, b0, sh, sb, aL, aR, A, offset, vn, ds)
    update_vs_flux!(fh, fb, bit_L, vs_data, vs_data_n, offset, rot)
end
function calc_flux!(
    ::Val{2},
    ::Val{-1},
    face::Face{T},
    DVM_data::DVM_Data,
    faceid::Int,
) where {T<:Ghost_PS_Data}
    ps_data = face.data
    dir = get_dir(faceid)
    rot = get_rot(faceid) # Left: -1, Right: 1
    gas = DVM_data.global_data.gas
    ds = ps_data.ds[dir]
    vs_data = ps_data.vs_data
    nps_data = ps_data.neighbor.data[faceid][1]
    vs_data_n = nps_data.vs_data
    weight_n, midpoint_n, df_n, sdf_n, offset, bit_L, _ =
        make_face_data(ps_data, nps_data, rot, dir)
    h = @view(df_n[:, 1])
    b = @view(df_n[:, 2])
    u = @view(midpoint_n[:, 1])
    v = @view(midpoint_n[:, 2])
    vn = @view(midpoint_n[:, dir])
    sh = @view(sdf_n[:, 1])
    sb = @view(sdf_n[:, 2])
    h0, b0 = reconstruct_vs(h, b, sh, sb, ds, 2.0 * ds, offset, rot)
    w0 = calc_w0(h0, b0, u, v, weight_n)
    prim0 = get_prim(w0, gas.γ)
    qf0 = calc_qf(h0, b0, prim0, u, v, weight_n)
    aL, aR = calc_a(w0, prim0, ps_data.w, nps_data.w, ds, 2.0 * ds, gas.K, rot)
    Mu, Mv, Mξ, Mu_L, Mu_R = moment_u(prim0, 6, 4, gas.K, rot, dir)
    A = calc_A(prim0, aL, aR, Mu, Mv, Mu_L, Mu_R, Mξ, gas.K, dir)
    τ0 = get_τ(prim0, gas.μᵣ, gas.ω)
    Mt = calc_time_int(τ0, gas.Δt)
    fw = calc_flux_g0(prim0, Mt, Mu, Mu_L, Mu_R, Mv, Mξ, aL, aR, A, dir)
    H, B = discrete_maxwell(u, v, prim0, gas.K)
    H⁺, B⁺ = shakhov_part(H, B, prim0, u, v, qf0, gas.Pr, gas.K)
    fw .+= calc_flux_f0(Mt, H⁺, B⁺, u, v, weight_n, h0, b0, sh, sb, vn)
    fw .*= ds * rot
    ps_data.flux .+= fw
    fh, fb = calc_fhb(H, B, H⁺, B⁺, Mt, Mξ, u, v, h0, b0, sh, sb, aL, aR, A, offset, vn, ds)
    update_vs_flux!(fh, fb, bit_L, vs_data, vs_data_n, offset, rot)
end

function update_flux!(DVM_data::DVM_Data)
    faces = DVM_data.faces
    @inbounds @simd for i in eachindex(faces)
        face = faces[i]
        faceid = face.faceid
        calc_flux!(
            Val(face.bound),
            Val(face.data.neighbor.state[faceid]),
            face,
            DVM_data,
            faceid,
        )
    end
end
