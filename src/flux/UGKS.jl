function calc_domain_flux(::UGKS,here_vs::Face_VS_Data,face::DomainFace{2,2,SuperSonicInflow},amr::AMR)
    rot,direction,midpoint,_,here_data = unpack(face)
    heavi,_,_,_,here_df,here_sdf = unpack(here_vs)
    global_data = amr.global_data
    gas = global_data.config.gas
    Δt = global_data.status.Δt
    vs_data = here_data.vs_data
    nheavi = [!x for x in heavi]
    there_mid = @views vs_data.midpoint[nheavi,:]
    there_weight = @views vs_data.weight[nheavi]
    there_vn = @views there_mid[:,direction]
    there_sdf = @views vs_data.sdf[nheavi,:,:]
    dx = midpoint-here_data.midpoint
    bc = get_bc(face.domain.bc)
    domain_w = get_conserved(bc,global_data)
    @inbounds @views begin
        df = [here_df[i,j]+dot(dx,here_sdf[i,j,:]) for i in axes(here_df,1),j in axes(here_df,2)]
    end
    ndf = Matrix{Float64}(undef,size(there_mid,1),2)
    @inbounds for i in axes(ndf,1)
        @views ndf[i,:] .= discrete_maxwell(there_mid[i,:],bc,global_data)
    end
    Here_vs = Face_VS_Data(here_vs,df);There_vs = Face_VS_Data{2,2}(nheavi,there_weight,there_mid,there_vn,ndf,there_sdf)
    w0 = calc_w0(Here_vs,There_vs)
    prim0 = get_prim(w0, global_data)
    qf0 = calc_qf(Here_vs,There_vs,prim0)
    aL, aR = calc_a(w0, prim0, here_data.w, domain_w, here_data.ds[direction], here_data.ds[direction], global_data, rot)
    Mu, Mv, Mξ, Mu_L, Mu_R = moment_u(prim0, global_data, rot, direction)
    A = calc_A(prim0, aL, aR, Mu, Mv, Mξ, Mu_L, Mu_R, global_data, direction)
    τ0 = get_τ(prim0, gas.μᵣ, gas.ω)
    Mt = time_int(τ0, Δt)
    fw = calc_flux_g0_2D2F(prim0, Mt, Mu, Mv, Mξ, Mu_L, Mu_R, aL, aR, A, direction)
    here_F = discrete_maxwell(Here_vs.midpoint, prim0, global_data)
    there_F = discrete_maxwell(There_vs.midpoint, prim0, global_data)
    here_F⁺ = shakhov_part(Here_vs.midpoint, here_F, prim0, qf0, global_data)
    there_F⁺ = shakhov_part(There_vs.midpoint, there_F, prim0, qf0, global_data)
    fw += calc_flux_f0(Here_vs, There_vs, here_F⁺, there_F⁺, Mt, direction)
    here_micro = calc_micro_flux(Here_vs, here_F, here_F⁺, aL, A, Mξ, Mt, direction)
    there_micro = calc_micro_flux(There_vs, there_F, there_F⁺, aR, A, Mξ, Mt, direction)
    return fw,[here_micro,there_micro]
end
function calc_domain_flux(::UGKS,here_vs::Face_VS_Data,face::DomainFace{2,2,UniformOutflow},amr::AMR)
    rot,direction,midpoint,_,here_data = unpack(face)
    heavi,_,_,_,here_df,here_sdf = unpack(here_vs)
    global_data = amr.global_data
    gas = global_data.config.gas
    Δt = global_data.status.Δt
    vs_data = here_data.vs_data
    nheavi = [!x for x in heavi]
    there_mid = @views vs_data.midpoint[nheavi,:]
    there_weight = @views vs_data.weight[nheavi]
    there_vn = @views there_mid[:,direction]
    there_df = @views vs_data.df[nheavi,:]
    there_sdf = @views vs_data.sdf[nheavi,:,:]
    dx = midpoint-here_data.midpoint
    @inbounds @views begin
        df = [here_df[i,j]+dot(dx,here_sdf[i,j,:]) for i in axes(here_df,1),j in axes(here_df,2)]
        ndf = [there_df[i,j]+dot(dx,there_sdf[i,j,:]) for i in axes(there_df,1),j in axes(there_df,2)]
    end
    Here_vs = Face_VS_Data(here_vs,df);There_vs = Face_VS_Data{2,2}(nheavi,there_weight,there_mid,there_vn,ndf,zeros(size(there_sdf)))
    w0 = calc_w0(Here_vs,There_vs)
    prim0 = get_prim(w0, global_data)
    qf0 = calc_qf(Here_vs,There_vs,prim0)
    aL, aR = calc_a(w0, prim0, here_data.w, here_data.w, here_data.ds[direction], here_data.ds[direction], global_data, rot)
    Mu, Mv, Mξ, Mu_L, Mu_R = moment_u(prim0, global_data, rot, direction)
    A = calc_A(prim0, aL, aR, Mu, Mv, Mξ, Mu_L, Mu_R, global_data, direction)
    τ0 = get_τ(prim0, gas.μᵣ, gas.ω)
    Mt = time_int(τ0, Δt)
    fw = calc_flux_g0_2D2F(prim0, Mt, Mu, Mv, Mξ, Mu_L, Mu_R, aL, aR, A, direction)
    here_F = discrete_maxwell(Here_vs.midpoint, prim0, global_data)
    there_F = discrete_maxwell(There_vs.midpoint, prim0, global_data)
    here_F⁺ = shakhov_part(Here_vs.midpoint, here_F, prim0, qf0, global_data)
    there_F⁺ = shakhov_part(There_vs.midpoint, there_F, prim0, qf0, global_data)
    fw += calc_flux_f0(Here_vs, There_vs, here_F⁺, there_F⁺, Mt, direction)
    here_micro = calc_micro_flux(Here_vs, here_F, here_F⁺, aL, A, Mξ, Mt, direction)
    there_micro = calc_micro_flux(There_vs, there_F, there_F⁺, aR, A, Mξ, Mt, direction)
    return fw,[here_micro,there_micro]
end
function calc_w0(here_vs::Face_VS_Data{2,2},there_vs::Face_VS_Data)
    @inbounds @views micro_to_macro_2D2F(
        here_vs.midpoint[:,1],here_vs.midpoint[:,2],
        here_vs.df[:,1],here_vs.df[:,2],here_vs.weight
    ) +
    micro_to_macro_2D2F(
        there_vs.midpoint[:,1],there_vs.midpoint[:,2],
        there_vs.df[:,1],there_vs.df[:,2],there_vs.weight
    )
end
function calc_qf(here_vs::Face_VS_Data{2,2},there_vs::Face_VS_Data,prim::AbstractVector)
    @inbounds @views heat_flux_2D2F(
        here_vs.midpoint[:,1],here_vs.midpoint[:,2],
        here_vs.df[:,1],here_vs.df[:,2],
        prim,here_vs.weight,
    ) +
    heat_flux_2D2F(
        there_vs.midpoint[:,1],there_vs.midpoint[:,2],
        there_vs.df[:,1],there_vs.df[:,2],
        prim,there_vs.weight,
    )
end
function calc_flux_f0(hvd::Face_VS_Data{2,2},tvd,here_F::AbstractMatrix,there_F,Mt,direction::Int)
    @inbounds @views calc_flux_f0_2D2F(
        hvd.midpoint[:,1],hvd.midpoint[:,2],
        hvd.df[:,1],hvd.df[:,2],
        hvd.sdf[:,1,direction],hvd.sdf[:,2,direction],
        hvd.weight,here_F[:,1],here_F[:,2],hvd.vn,Mt
    ) +
    calc_flux_f0_2D2F(
        tvd.midpoint[:,1],tvd.midpoint[:,2],
        tvd.df[:,1],tvd.df[:,2],
        tvd.sdf[:,1,direction],tvd.sdf[:,2,direction],
        tvd.weight,there_F[:,1],there_F[:,2],tvd.vn,Mt
    )
end
function calc_micro_flux_2D2F(
    u::AbstractVector{T},
    v::AbstractVector,
    h::AbstractVector,
    b::AbstractVector,
    sh::AbstractVector,
    sb::AbstractVector,
    H0::AbstractVector,
    B0::AbstractVector,
    H⁺::AbstractVector,
    B⁺::AbstractVector,
    ax::AbstractVector,
    at::AbstractVector,
    Mξ::AbstractVector,
    Mt::AbstractVector,
    vn::AbstractVector,
) where {T}
    micro_flux = Matrix{T}(undef, length(u), 2)
    @inbounds begin
        @. @view(micro_flux[:, 1]) =
            Mt[1] * vn * (H0 + H⁺) +
            Mt[2] *
            vn^2 *
            (
                ax[1] * H0 +
                ax[2] * u * H0 +
                ax[3] * v * H0 +
                0.5 * ax[4] * ((u^2 + v^2) * H0 + B0)
            )+
            Mt[3] *
            vn *
            (
                at[1] * H0 +
                at[2] * u * H0 +
                at[3] * v * H0 +
                0.5 * at[4] * ((u^2 + v^2) * H0 + B0)
            ) +
            Mt[4] * vn * h - Mt[5] * vn^2 * sh
        @. @view(micro_flux[:, 2]) =
            Mt[1] * vn * (B0 + B⁺) +
            Mt[2] *
            vn^2 *
            (
                ax[1] * B0 +
                ax[2] * u * B0 +
                ax[3] * v * B0 +
                0.5 * ax[4] * ((u^2 + v^2) * B0 + Mξ[3] * H0)
            ) +
            Mt[3] *
            vn *
            (
                at[1] * B0 +
                at[2] * u * B0 +
                at[3] * v * B0 +
                0.5 * at[4] * ((u^2 + v^2) * B0 + Mξ[3] * H0)
            ) +
            Mt[4] * vn * b - Mt[5] * vn^2 * sb
    end
    return micro_flux
end
function calc_micro_flux(
    f_vs_data::Face_VS_Data{2,2},
    F0::AbstractMatrix{T},
    F⁺,
    ax,
    at,
    Mξ,
    Mt,
    direction,
) where {T}
    @inbounds @views begin
        calc_micro_flux_2D2F(
            f_vs_data.midpoint[:, 1],
            f_vs_data.midpoint[:, 2],
            f_vs_data.df[:, 1],
            f_vs_data.df[:, 2],
            f_vs_data.sdf[:, 1,direction],
            f_vs_data.sdf[:, 2,direction],
            F0[:, 1],
            F0[:, 2],
            F⁺[:, 1],
            F⁺[:, 2],
            ax,
            at,
            Mξ,
            Mt,
            f_vs_data.vn,
        )
    end
end
function calc_flux(::UGKS,here_vs,there_vs,flux_data::Union{FullFace,Flux_Data},amr::AMR{2,2}) # without face area
    rot,direction,midpoint,here_data,there_data = unpack(flux_data)
    global_data = amr.global_data
    gas = global_data.config.gas;Δt = global_data.status.Δt
    here_df = here_vs.df;there_df = there_vs.df
    here_sdf = here_vs.sdf;there_sdf = there_vs.sdf
    here_ps_mid = here_data.midpoint;there_ps_mid = there_data.midpoint
    dx = midpoint-here_ps_mid
    ndx = midpoint-there_ps_mid
    @inbounds @views begin
        df = [here_df[i,j]+dot(dx,here_sdf[i,j,:]) for i in axes(here_df,1),j in axes(here_df,2)]
        ndf = [there_df[i,j]+dot(ndx,there_sdf[i,j,:]) for i in axes(there_df,1),j in axes(there_df,2)]
    end
    Here_vs = Face_VS_Data(here_vs,df);There_vs = Face_VS_Data(there_vs,ndf)
    w0 = calc_w0(Here_vs,There_vs)
    prim0 = get_prim(w0, global_data)
    qf0 = calc_qf(Here_vs,There_vs,prim0)
    wL = @. here_data.w+dx[FAT[1][direction]]*here_data.sw[:,FAT[1][direction]]
    wR = @. there_data.w+ndx[FAT[1][direction]]*there_data.sw[:,FAT[1][direction]]
    aL, aR = calc_a(w0, prim0, wL, wR, here_data.ds[direction], there_data.ds[direction], global_data, rot)
    Mu, Mv, Mξ, Mu_L, Mu_R = moment_u(prim0, global_data, rot, direction)
    A = calc_A(prim0, aL, aR, Mu, Mv, Mξ, Mu_L, Mu_R, global_data, direction)
    τ0 = get_τ(prim0, gas.μᵣ, gas.ω)
    Mt = time_int(τ0, Δt)
    fw = calc_flux_g0_2D2F(prim0, Mt, Mu, Mv, Mξ, Mu_L, Mu_R, aL, aR, A, direction)
    here_F = discrete_maxwell(Here_vs.midpoint, prim0, global_data)
    there_F = discrete_maxwell(There_vs.midpoint, prim0, global_data)
    here_F⁺ = shakhov_part(Here_vs.midpoint, here_F, prim0, qf0, global_data)
    there_F⁺ = shakhov_part(There_vs.midpoint, there_F, prim0, qf0, global_data)
    fw += calc_flux_f0(Here_vs, There_vs, here_F⁺, there_F⁺, Mt, direction)
    here_micro = calc_micro_flux(Here_vs, here_F, here_F⁺, aL, A, Mξ, Mt, direction)
    there_micro = calc_micro_flux(There_vs, there_F, there_F⁺, aR, A, Mξ, Mt, direction)
    return fw,[here_micro,there_micro]
end