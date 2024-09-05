get_dir(faceid::Int) = div(faceid - 1, 2) + 1
get_rot(faceid::Int) = (-1.0)^(faceid - 1)
get_sound(prim::Vector{Cdouble}, γ::Real) = √(0.5 * γ * prim[end])
heaviside(x::Real) = ifelse(x >= 0, one(x), zero(x))
function get_prim(ps_data::PS_Data{2,NDF}, global_data::Global_Data{2,NDF}) where {NDF}
    get_prim_2D(ps_data.w, global_data.config.gas.γ)
end
function get_prim(w::AbstractVector, global_data::Global_Data{2,NDF}) where {NDF}
    get_prim_2D(w, global_data.config.gas.γ)
end
function get_prim(ps_data::PS_Data{3,NDF}, global_data::Global_Data{3,NDF}) where {NDF}
    get_prim_3D(ps_data.w, global_data.config.gas.γ)
end
function get_prim(w::AbstractVector, global_data::Global_Data{3,NDF}) where {NDF}
    get_prim_3D(w, global_data.config.gas.γ)
end

function get_conserved(ps_data::PS_Data{3,NDF}, global_data::Global_Data) where {NDF}
    get_conserved_3D(ps_data.prim, global_data.config.gas.γ)
end

function discrete_maxwell(ps_data::PS_Data{3,1}, global_data::Global_Data)
    discrete_maxwell(ps_data.vs_data.midpoint, ps_data.prim, global_data)
end
function discrete_maxwell(
    midpoint::AbstractMatrix,
    prim::AbstractVector,
    ::Global_Data{3,1},
)
    discrete_maxwell_3D1F(
        @view(midpoint[:, 1]),
        @view(midpoint[:, 2]),
        @view(midpoint[:, 3]),
        prim,
    )
end
function discrete_maxwell(ps_data::PS_Data{2,2}, global_data::Global_Data)
    discrete_maxwell(ps_data.vs_data.midpoint, ps_data.prim, global_data)
end
function discrete_maxwell(
    midpoint::AbstractMatrix,
    prim::AbstractVector,
    ::Global_Data{2,2},
)
    discrete_maxwell_2D2F(
        @view(midpoint[:, 1]),
        @view(midpoint[:, 2]),
        prim,
        global_data.config.gas.K,
    )
end
function shakhov_part(
    midpoint::AbstractMatrix,
    F::AbstractMatrix,
    prim::AbstractVector,
    qf::AbstractVector,
    global_data::Global_Data{2,2},
)
    shakhov_part_2D2F(
        @view(midpoint[:, 1]),
        @view(midpoint[:, 2]),
        @view(F[:, 1]),
        @view(F[:, 2]),
        prim,
        qf,
        global_data.config.gas.Pr,
        global_data.config.gas.K,
    )
end
function shakhov_part(
    midpoint::AbstractMatrix,
    F::AbstractMatrix,
    prim::AbstractVector,
    qf::AbstractVector,
    global_data::Global_Data{3,1},
)
    shakhov_part_3D1F(
        @view(midpoint[:, 1]),
        @view(midpoint[:, 2]),
        @view(midpoint[:, 3]),
        vec(F),
        prim,
        qf,
        global_data.config.gas.Pr,
    )
end
function shakhov_part(
    midpoint::AbstractMatrix,
    F::AbstractVector,
    prim::AbstractVector,
    qf::AbstractVector,
    global_data::Global_Data{3,1},
)
    shakhov_part_3D1F(
        @view(midpoint[:, 1]),
        @view(midpoint[:, 2]),
        @view(midpoint[:, 3]),
        F,
        prim,
        qf,
        global_data.config.gas.Pr,
    )
end

function calc_qf(vs_data::AbstractVsData{2,2}, prim::AbstractVector)
    heat_flux_2D2F(
        @view(vs_data.midpoint[:, 1]),
        @view(vs_data.midpoint[:, 2]),
        @view(vs_data.df[:, 1]),
        @view(vs_data.df[:, 2]),
        prim,
        vs_data.weight,
    )
end
function calc_qf(vs_data::AbstractVsData{3,1}, prim::AbstractVector)
    heat_flux_3D1F(
        @view(vs_data.midpoint[:, 1]),
        @view(vs_data.midpoint[:, 2]),
        @view(vs_data.midpoint[:, 3]),
        vec(vs_data.df),
        prim,
        vs_data.weight,
    )
end

function calc_ρw(
    vs_data::VS_Data{2,2},
    df0::AbstractMatrix,
    prim0::AbstractVector,
    Θ::AbstractVector,
    vn::AbstractVector,
)
    maxwellian_density_2D2F(
        @view(vs_data.midpoint[:, 1]),
        @view(vs_data.midpoint[:, 2]),
        @view(df0[:, 1]),
        prim0,
        vs_data.weight,
        Θ,
        vn,
    )
end
function calc_ρw(
    vs_data::VS_Data{3,1},
    f0::AbstractMatrix,
    prim0::AbstractVector,
    Θ::AbstractVector,
    vn::AbstractVector,
)
    maxwellian_density_3D1F(
        @view(vs_data.midpoint[:, 1]),
        @view(vs_data.midpoint[:, 2]),
        @view(vs_data.midpoint[:, 3]),
        vec(f0),
        prim0,
        vs_data.weight,
        Θ,
        vn,
    )
end

function calc_fwb(vs_data::VS_Data{2,2}, F::AbstractMatrix, vn::AbstractVector)
    @inbounds macro_flux_2D2F(
        @view(vs_data.midpoint[:, 1]),
        @view(vs_data.midpoint[:, 2]),
        @view(F[:, 1]),
        @view(F[:, 2]),
        vs_data.weight,
        vn,
    )
end
function calc_fwb(vs_data::VS_Data{3,1}, F::AbstractMatrix, vn::AbstractVector)
    @inbounds macro_flux_3D1F(
        @view(vs_data.midpoint[:, 1]),
        @view(vs_data.midpoint[:, 2]),
        @view(vs_data.midpoint[:, 3]),
        vec(F),
        vs_data.weight,
        vn,
    )
end

function reconstruct_vs!(
    f_vs_data::Face_VS_Data,
    dsL::T,
    dsR::Real,
    offset::Integer,
    rot::Real,
) where {T}
    df_n = f_vs_data.df
    sdf_n = f_vs_data.sdf
    @inbounds @. df_n[1:offset, :] -= @views 0.5 * rot * dsL * sdf_n[1:offset, :]
    @inbounds @. df_n[offset+1:end, :] += @views 0.5 * rot * dsR * sdf_n[offset+1:end, :]
end
function reconstruct_vs_L!(f_vs_data::Face_VS_Data, dsL::Real, offset::Integer, rot::Real)
    @inbounds @. f_vs_data.df[1:offset, :] -=
        @views 0.5 * rot * dsL * f_vs_data.sdf[1:offset, :]
end
function reconstruct_vs_R!(f_vs_data::Face_VS_Data, dsR::Real, offset::Integer, rot::Real)
    @inbounds @. f_vs_data.df[offset+1:end, :] +=
        @views 0.5 * rot * dsR * f_vs_data.sdf[offset+1:end, :]
end

function calc_w0(f_vs_data::Face_VS_Data{2,2})
    @inbounds micro_to_macro_2D2F(
        @view(f_vs_data.midpoint[:, 1]),
        @view(f_vs_data.midpoint[:, 2]),
        @view(f_vs_data.df[:, 1]),
        @view(f_vs_data.df[:, 2]),
        f_vs_data.weight,
    )
end
function calc_w0(f_vs_data::Face_VS_Data{3,1})
    @inbounds micro_to_macro_3D1F(
        @view(f_vs_data.midpoint[:, 1]),
        @view(f_vs_data.midpoint[:, 2]),
        @view(f_vs_data.midpoint[:, 3]),
        vec(f_vs_data.df[:, 1]),
        f_vs_data.weight,
    )
end
function calc_a(
    w0::AbstractVector{T},
    prim0::AbstractVector,
    cLw::AbstractVector,
    cRw,
    dxL::Real,
    dxR,
    global_data::Global_Data,
    rot::Real,
) where {T}
    sw = Vector{T}(undef, 5)
    @. sw = rot * (cLw - w0) / (0.5 * dxL)
    aL = micro_slope(sw, prim0, global_data)
    @. sw = rot * (w0 - cRw) / (0.5 * dxR)
    aR = micro_slope(sw, prim0, global_data)
    (aL, aR)
end

function micro_slope(sw::AbstractVector, prim::AbstractVector, global_data::Global_Data{2})
    micro_slope_2D(sw, prim, global_data.config.gas.K)
end
function micro_slope(sw::AbstractVector, prim::AbstractVector, global_data::Global_Data{3})
    micro_slope_3D(sw, prim, global_data.config.gas.K)
end

function moment_u(
    prim0::AbstractVector,
    global_data::Global_Data{2,2},
    ROT::Real,
    DIR::Integer,
)
    @inbounds U = prim0[DIR+1]
    V = prim0[FAT[1][DIR]+1]
    λ = prim0[4]
    Mu, Mv, Mξ, Mu_L, Mu_R = moment_u_2D2F(U, V, λ, 6, 4, global_data.config.gas.K)
    ROT < 0 && return (Mu, Mv, Mξ, Mu_L, Mu_R)
    return (Mu, Mv, Mξ, Mu_R, Mu_L)
end
function moment_u(prim0::AbstractVector, ::Global_Data{3,1}, ROT::Real, DIR::Integer)
    @inbounds begin 
        U = prim0[DIR+1]
        V = prim0[FAT[2][DIR][1]+1]
        W = prim0[FAT[2][DIR][2]+1]
        λ = prim0[5]
    end
    Mu, Mv, Mw, Mu_L, Mu_R = moment_u_3D1F(U, V, W, λ, 6, 4, 4)
    ROT < 0 && return (Mu, Mv, Mw, Mu_L, Mu_R)
    return (Mu, Mv, Mw, Mu_R, Mu_L)
end

function calc_A(
    prim::AbstractVector,
    aL,
    aR,
    Mu,
    Mv,
    Mξ,
    Mu_L,
    Mu_R,
    global_data::Global_Data{2,2},
    dir::Int,
)
    if dir == 1
        Mau_L = moment_au_2D2F(aL, Mu_L, Mv, Mξ, 1, 0)
        Mau_R = moment_au_2D2F(aR, Mu_R, Mv, Mξ, 1, 0)
    else
        Mau_L = moment_au_2D2F(aL, Mu, Mu_L, Mξ, 0, 1)
        Mau_R = moment_au_2D2F(aR, Mu, Mu_R, Mξ, 0, 1)
    end
    @inbounds sw = -prim[1] * (Mau_L + Mau_R)
    micro_slope(sw, prim, global_data)
end
function calc_A(
    prim::AbstractVector,
    aL,
    aR,
    Mu,
    Mv,
    Mw,
    Mu_L,
    Mu_R,
    global_data::Global_Data{3,1},
    dir::Int,
)
    if dir == 1
        Mau_L = moment_au_3D1F(aL, Mu_L, Mv, Mw, 1, 0, 0)
        Mau_R = moment_au_3D1F(aR, Mu_R, Mv, Mw, 1, 0, 0)
    elseif dir == 2
        Mau_L = moment_au_3D1F(aL, Mu, Mu_L, Mw, 0, 1, 0)
        Mau_R = moment_au_3D1F(aR, Mu, Mu_R, Mw, 0, 1, 0)
    else
        Mau_L = moment_au_3D1F(aL, Mu, Mv, Mu_L, 0, 0, 1)
        Mau_R = moment_au_3D1F(aR, Mu, Mv, Mu_R, 0, 0, 1)
    end
    @inbounds sw = -prim[1] * (Mau_L + Mau_R)
    micro_slope(sw, prim, global_data)
end

function calc_flux_g0_2D2F(prim, Mt, Mu, Mv, Mξ, Mu_L, Mu_R, aL, aR, A, dir)
    if dir == 1
        Mau_0 = moment_uv_2D2F(Mu, Mv, Mξ, 1, 0, 0)
        Mau2_L = moment_au_2D2F(aL, Mu_L, Mv, Mξ, 2, 0)
        Mau2_R = moment_au_2D2F(aR, Mu_R, Mv, Mξ, 2, 0)
        Mau_T = moment_au_2D2F(A, Mu, Mv, Mξ, 1, 0)
    else
        Mau_0 = moment_uv_2D2F(Mu, Mv, Mξ, 0, 1, 0)
        Mau2_L = moment_au_2D2F(aL, Mu, Mu_L, Mξ, 0, 2)
        Mau2_R = moment_au_2D2F(aR, Mu, Mu_R, Mξ, 0, 2)
        Mau_T = moment_au_2D2F(A, Mu, Mv, Mξ, 0, 1)
    end
    @inbounds prim[1] * (Mt[1] * Mau_0 + Mt[2] * (Mau2_L + Mau2_R) + Mt[3] * Mau_T)
end

function calc_flux_g0_3D1F(prim, Mt, Mu, Mv, Mw, Mu_L, Mu_R, aL, aR, A, dir)
    if dir == 1
        Mau_0 = moment_uv_3D1F(Mu, Mv, Mw, 1, 0, 0)
        Mau2_L = moment_au_3D1F(aL, Mu_L, Mv, Mw, 2, 0, 0)
        Mau2_R = moment_au_3D1F(aR, Mu_R, Mv, Mw, 2, 0, 0)
        Mau_T = moment_au_3D1F(A, Mu, Mv, Mw, 1, 0, 0)
    elseif dir == 2
        Mau_0 = moment_uv_3D1F(Mu, Mv, Mw, 0, 1, 0)
        Mau2_L = moment_au_3D1F(aL, Mu, Mu_L, Mw, 0, 2, 0)
        Mau2_R = moment_au_3D1F(aR, Mu, Mu_R, Mw, 0, 2, 0)
        Mau_T = moment_au_3D1F(A, Mu, Mv, Mw, 0, 1, 0)
    else
        Mau_0 = moment_uv_3D1F(Mu, Mv, Mw, 0, 0, 1)
        Mau2_L = moment_au_3D1F(aL, Mu, Mu_L, Mw, 0, 0, 2)
        Mau2_R = moment_au_3D1F(aR, Mu, Mu_R, Mw, 0, 0, 2)
        Mau_T = moment_au_3D1F(A, Mu, Mv, Mw, 0, 0, 1)
    end
    @inbounds prim[1] * (Mt[1] * Mau_0 + Mt[2] * (Mau2_L + Mau2_R) + Mt[3] * Mau_T)
end

function calc_flux_f0(f_vs_data::Face_VS_Data{2,2}, F⁺::AbstractMatrix, Mt::AbstractVector)
    @inbounds calc_flux_f0_2D2F(
        @view(f_vs_data.midpoint[:, 1]),
        @view(f_vs_data.midpoint[:, 2]),
        @view(f_vs_data.df[:, 1]),
        @view(f_vs_data.df[:, 2]),
        @view(f_vs_data.sdf[:, 1]),
        @view(f_vs_data.sdf[:, 2]),
        f_vs_data.weight,
        @view(F⁺[:, 1]),
        @view(F⁺[:, 2]),
        f_vs_data.vn,
        Mt,
    )
end
function calc_flux_f0(f_vs_data::Face_VS_Data{3,1}, F⁺::AbstractVector, Mt::AbstractVector)
    @inbounds calc_flux_f0_3D1F(
        @view(f_vs_data.midpoint[:, 1]),
        @view(f_vs_data.midpoint[:, 2]),
        @view(f_vs_data.midpoint[:, 3]),
        vec(f_vs_data.df),
        f_vs_data.sdf,
        f_vs_data.weight,
        F⁺,
        f_vs_data.vn,
        Mt,
    )
end

function calc_flux_f0_2D2F(u::AbstractVector{T}, v, h, b, sh, sb, weight, H⁺, B⁺, vn, Mt) where {T}
    F = Vector{T}(undef, 4)
    @inbounds begin
        F[1] =
            Mt[1] * sum(weight .* vn .* H⁺) + Mt[4] * sum(weight .* vn .* h) -
            Mt[5] * sum(weight .* vn .^ 2 .* sh)
        F[2] =
            Mt[1] * sum(weight .* vn .* u .* H⁺) + Mt[4] * sum(weight .* vn .* u .* h) -
            Mt[5] * sum(weight .* vn .^ 2 .* u .* sh)
        F[3] =
            Mt[1] * sum(weight .* vn .* v .* H⁺) + Mt[4] * sum(weight .* vn .* v .* h) -
            Mt[5] * sum(weight .* v .* vn .^ 2 .* sh)
        F[4] =
            Mt[1] *
            0.5 *
            (sum(weight .* vn .* (u .^ 2 + v .^ 2) .* H⁺) + sum(weight .* vn .* B⁺)) +
            Mt[4] *
            0.5 *
            (sum(weight .* vn .* (u .^ 2 + v .^ 2) .* h) + sum(weight .* vn .* b)) -
            Mt[5] *
            0.5 *
            (
                sum(weight .* vn .^ 2 .* (u .^ 2 + v .^ 2) .* sh) +
                sum(weight .* vn .^ 2 .* sb)
            )
    end
    return F
end

function calc_flux_f0_3D1F(u::AbstractVector{T}, v, w, f, sf, weight, F⁺, vn, Mt) where {T}
    F = Vector{T}(undef, 5)
    @inbounds begin
        F[1] =
            Mt[1] * sum(weight .* vn .* F⁺) + Mt[4] * sum(weight .* vn .* f) -
            Mt[5] * sum(weight .* vn .^ 2 .* sf)
        F[2] =
            Mt[1] * sum(weight .* vn .* u .* F⁺) + Mt[4] * sum(weight .* vn .* u .* f) -
            Mt[5] * sum(weight .* vn .^ 2 .* u .* sf)
        F[3] =
            Mt[1] * sum(weight .* vn .* v .* F⁺) + Mt[4] * sum(weight .* vn .* v .* f) -
            Mt[5] * sum(weight .* v .* vn .^ 2 .* sf)
        F[4] =
            Mt[1] * sum(weight .* vn .* w .* F⁺) + Mt[4] * sum(weight .* vn .* w .* f) -
            Mt[5] * sum(weight .* w .* vn .^ 2 .* sf)
        F[5] =
            Mt[1] * 0.5 * sum(weight .* vn .* (u .^ 2 + v .^ 2 + w .^ 2) .* F⁺) +
            Mt[4] * 0.5 * sum(weight .* vn .* (u .^ 2 + v .^ 2 + w .^ 2) .* f) -
            Mt[5] * 0.5 * sum(weight .* vn .^ 2 .* (u .^ 2 + v .^ 2 + w .^ 2) .* sf)
    end
    return F
end

function calc_micro_flux(
    f_vs_data::Face_VS_Data{2,2},
    F0::AbstractMatrix{T},
    F⁺,
    aL,
    aR,
    aT,
    Mξ,
    Mt,
    offset,
    dsf,
) where {T}
    Θ = Array{T}(undef, size(F0, 1))
    @inbounds begin
        Θ[1:offset] .= 1.0
        Θ[offset+1:end] .= 0.0
        calc_micro_flux_2D2F(
            @view(f_vs_data.midpoint[:, 1]),
            @view(f_vs_data.midpoint[:, 2]),
            @view(f_vs_data.df[:, 1]),
            @view(f_vs_data.df[:, 2]),
            @view(f_vs_data.sdf[:, 1]),
            @view(f_vs_data.sdf[:, 2]),
            @view(F0[:, 1]),
            @view(F0[:, 2]),
            @view(F⁺[:, 1]),
            @view(F⁺[:, 2]),
            aL,
            aR,
            aT,
            Mξ,
            Mt,
            f_vs_data.vn,
            Θ,
            dsf,
        )
    end
end
function calc_micro_flux(
    f_vs_data::Face_VS_Data{3,1},
    F0::AbstractVector{T},
    F⁺,
    aL,
    aR,
    aT,
    Mt,
    offset,
    dsf,
) where {T}
    Θ = Array{T}(undef, size(F0, 1))
    @inbounds begin
        Θ[1:offset] .= 1.0
        Θ[offset+1:end] .= 0.0
        calc_micro_flux_3D1F(
            @view(f_vs_data.midpoint[:, 1]),
            @view(f_vs_data.midpoint[:, 2]),
            @view(f_vs_data.midpoint[:, 3]),
            vec(f_vs_data.df),
            vec(f_vs_data.sdf),
            F0,
            F⁺,
            aL,
            aR,
            aT,
            Mt,
            f_vs_data.vn,
            Θ,
            dsf,
        )
    end
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
    aL::AbstractVector,
    aR::AbstractVector,
    aT::AbstractVector,
    Mξ::AbstractVector,
    Mt::AbstractVector,
    vn::AbstractVector,
    Θ::AbstractVector,
    ds::Real,
) where {T}
    micro_flux = Matrix{T}(undef, length(u), 2)
    @inbounds begin
        @. @view(micro_flux[:, 1]) =
            Mt[1] * vn * (H0 + H⁺) +
            Mt[2] *
            vn^2 *
            (
                aL[1] * H0 +
                aL[2] * u * H0 +
                aL[3] * v * H0 +
                0.5 * aL[4] * ((u^2 + v^2) * H0 + B0)
            ) *
            Θ +
            Mt[2] *
            vn^2 *
            (
                aR[1] * H0 +
                aR[2] * u * H0 +
                aR[3] * v * H0 +
                0.5 * aR[4] * ((u^2 + v^2) * H0 + B0)
            ) *
            (-Θ + 1.0) +
            Mt[3] *
            vn *
            (
                aT[1] * H0 +
                aT[2] * u * H0 +
                aT[3] * v * H0 +
                0.5 * aT[4] * ((u^2 + v^2) * H0 + B0)
            ) +
            Mt[4] * vn * h - Mt[5] * vn^2 * sh
        @. @view(micro_flux[:, 2]) =
            Mt[1] * vn * (B0 + B⁺) +
            Mt[2] *
            vn^2 *
            (
                aL[1] * B0 +
                aL[2] * u * B0 +
                aL[3] * v * B0 +
                0.5 * aL[4] * ((u^2 + v^2) * B0 + Mξ[3] * H0)
            ) *
            Θ +
            Mt[2] *
            vn^2 *
            (
                aR[1] * B0 +
                aR[2] * u * B0 +
                aR[3] * v * B0 +
                0.5 * aR[4] * ((u^2 + v^2) * B0 + Mξ[3] * H0)
            ) *
            (-Θ + 1.0) +
            Mt[3] *
            vn *
            (
                aT[1] * B0 +
                aT[2] * u * B0 +
                aT[3] * v * B0 +
                0.5 * aT[4] * ((u^2 + v^2) * B0 + Mξ[3] * H0)
            ) +
            Mt[4] * vn * b - Mt[5] * vn^2 * sb
    end
    return micro_flux .* ds
end

function calc_micro_flux_3D1F(
    u::AbstractVector{T},
    v::AbstractVector,
    w::AbstractVector,
    f::AbstractVector,
    sf::AbstractVector,
    F0::AbstractVector,
    F⁺::AbstractVector,
    aL::AbstractVector,
    aR::AbstractVector,
    aT::AbstractVector,
    Mt::AbstractVector,
    vn::AbstractVector,
    Θ::AbstractVector,
    ds::Real,
) where {T}
    micro_flux = Matrix{T}(undef, length(u), 1)
    @inbounds begin
        @. micro_flux =
            Mt[1] * vn * (F0 + F⁺) +
            Mt[2] *
            vn^2 *
            (
                aL[1] * F0 +
                aL[2] * u * F0 +
                aL[3] * v * F0 +
                aL[4] * w * F0 +
                0.5 * aL[5] * (u^2 + v^2 + w^2) * F0
            ) *
            Θ +
            Mt[2] *
            vn^2 *
            (
                aR[1] * F0 +
                aR[2] * u * F0 +
                aR[3] * v * F0 +
                aR[4] * w * F0 +
                0.5 * aR[5] * (u^2 + v^2 + w^2) * F0
            ) *
            (-Θ + 1.0) +
            Mt[3] *
            vn *
            (
                aT[1] * F0 +
                aT[2] * u * F0 +
                aT[3] * v * F0 +
                aT[4] * w * F0 +
                0.5 * aT[5] * (u^2 + v^2) * F0
            ) +
            Mt[4] * vn * f - Mt[5] * vn^2 * sf
    end
    return micro_flux .* ds
end