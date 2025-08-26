get_dir(faceid::Int) = div(faceid - 1, 2) + 1
get_rot(faceid::Int) = (-1.0)^(faceid - 1)
get_sound(prim::Vector{Cdouble}, γ::Real) = √(0.5 * γ * prim[4])
heaviside(x::Real) = ifelse(x >= 0, one(x), zero(x))
function get_prim(w::Vector, γ::Real)
    prim = Vector{Float64}(undef, 4)
    prim[1] = w[1]
    prim[2] = w[2] / w[1]
    prim[3] = w[3] / w[1]
    prim[4] = 0.5 * w[1] / (γ - 1.0) / (w[4] - 0.5 * (w[2]^2 + w[3]^2) / w[1])
    prim
end
function get_conserved(prim::Vector, γ::Real)
    w = Vector{Float64}(undef, 4)
    w[1] = prim[1]
    w[2] = prim[1] * prim[2]
    w[3] = prim[1] * prim[3]
    w[4] = 0.5 * prim[1] / prim[4] / (γ - 1.0) + 0.5 * prim[1] * (prim[2]^2 + prim[3]^2)
    w
end
function localize(midpoint::AM, n::AV)
    n1 = [0.0 -1.0; 1.0 0.0]
    vn = midpoint * n
    vt = midpoint * n1 * n
    (vn, vt)
end
function localize(prim::AV, n::AV)
    p = Vector{Float64}(undef, 4)
    p[1] = prim[1]
    p[2] = prim[2] * n[1] + prim[3] * n[2]
    p[3] = -prim[2] * n[2] + prim[3] * n[1]
    p[4] = prim[4]
    p
end
function globalize(fw::AV, n::AV)
    f = Vector{Float64}(undef, 4)
    f[1] = fw[1]
    f[2] = fw[2] * n[1] - fw[3] * n[2]
    f[3] = fw[2] * n[2] + fw[3] * n[1]
    f[4] = fw[4]
    f
end
function calc_ρw(prim::AV, h::AV, Θ::AV, u::AV, v::AV, weight::AV, vn::AV)
    SF = sum(@. weight * vn * h * (1.0 - Θ))
    SG =
        prim[4] / π *
        sum(@. weight * vn * exp(-prim[4] * ((u - prim[2])^2 + (v - prim[3])^2)) * Θ)
    -SF / SG
end

function discrete_maxwell_2DF_2D(midpoint::AM, prim::AV, K::Real)
    DF = Array{Float64}(undef, size(midpoint, 1), 2)
    @. DF[:, 1] =
        prim[1] *
        (prim[4] / π) *
        exp.(
            -prim[4] *
            ((@view(midpoint[:, 1]) - prim[2])^2 + (@view(midpoint[:, 2]) - prim[3])^2),
        )
    @. DF[:, 2] = @view(DF[:, 1]) * K / (2.0 * prim[4])
    DF
end # 580.569 ns (9 allocations: 3.56 KiB) for vs_num=100
function discrete_maxwell_2DF_2D(u::AV, v::AV, prim::AV, K::Real)
    H = @. prim[1] * (prim[4] / π) * exp.(-prim[4] * ((u - prim[2])^2 + (v - prim[3])^2))
    B = @. H * K / (2.0 * prim[4])
    (H, B)
end
function shakhov_part_2DF_2D(DF::AM, prim::AV, midpoint::AM, qf::AV, Pr::Real, K::Real)
    u = @views midpoint[:, 1]
    v = @views midpoint[:, 2]
    H = @views DF[:, 1]
    B = @views DF[:, 2]
    H_plus = @. 0.8 * (1 - Pr) * prim[4]^2 / prim[1] *
       ((u - prim[2]) * qf[1] + (v - prim[3]) * qf[2]) *
       (2 * prim[4] * ((u - prim[2])^2 + (v - prim[3])^2) + K - 5) *
       H
    B_plus = @. 0.8 * (1 - Pr) * prim[4]^2 / prim[1] *
       ((u - prim[2]) * qf[1] + (v - prim[3]) * qf[2]) *
       (2 * prim[4] * ((u - prim[2])^2 + (v - prim[3])^2) + K - 3) *
       B
    (H_plus, B_plus)
end
function shakhov_part_2DF_2D(
    H::AV,
    B::AV,
    prim::AV,
    u::AV,
    v::AV,
    qf::AV,
    Pr::Real,
    K::Real,
)
    H_plus = @. 0.8 * (1 - Pr) * prim[4]^2 / prim[1] *
       ((u - prim[2]) * qf[1] + (v - prim[3]) * qf[2]) *
       (2 * prim[4] * ((u - prim[2])^2 + (v - prim[3])^2) + K - 5) *
       H
    B_plus = @. 0.8 * (1 - Pr) * prim[4]^2 / prim[1] *
       ((u - prim[2]) * qf[1] + (v - prim[3]) * qf[2]) *
       (2 * prim[4] * ((u - prim[2])^2 + (v - prim[3])^2) + K - 3) *
       B
    (H_plus, B_plus)
end
function calc_fwb_2DF_2D(df::AM, midpoint::AM, weight::AV, dir::Int)
    fw = zeros(DIM + 2)
    vn = @views midpoint[:, dir]
    h = @views df[:, 1]
    b = @views df[:, 2]
    u = @views midpoint[:, 1]
    v = @views midpoint[:, 2]
    fw[1] = sum(@. weight * vn * h)
    fw[2] = sum(@. weight * vn * u * h)
    fw[3] = sum(@. weight * vn * v * h)
    fw[4] = 0.5 * sum(@. weight * vn * ((u^2 + v^2) * h + b))
    fw
end
function calc_fwb_2DF_2D(h::AV, b::AV, u::AV, v::AV, weight::AV, vn::AV)
    fw = Array{Float64}(undef, 4)
    fw[1] = sum(@. weight * vn * h)
    fw[2] = sum(@. weight * vn * u * h)
    fw[3] = sum(@. weight * vn * v * h)
    fw[4] = 0.5 * sum(@. weight * vn * ((u^2 + v^2) * h + b))
    fw
end
function reconstruct_vs(
    df_n::AM,
    sdf_n::Matrix{Float64},
    dsL::Float64,
    dsR::Float64,
    offset::Int,
    rot::Float64,
)
    df0 = Array{Float64}(undef, size(df_n))
    @. df0[1:offset] = @views df_n[1:offset] - 0.5 * rot * dsL * sdf_n[1:offset]
    @. df0[(offset+1):end] =
        @views df_n[(offset+1):end] + 0.5 * rot * dsR * sdf_n[(offset+1):end]
    df0
end
function reconstruct_vs(
    h::AV,
    b,
    sh::AV,
    sb,
    dsL::Float64,
    dsR::Float64,
    offset::Int,
    rot::Float64,
)
    l = length(h)
    h0 = Vector{Float64}(undef, l)
    b0 = Vector{Float64}(undef, l)
    @. h0[1:offset] = @views h[1:offset] - 0.5 * rot * dsL * sh[1:offset]
    @. b0[1:offset] = @views b[1:offset] - 0.5 * rot * dsL * sb[1:offset]
    @. h0[(offset+1):end] = @views h[(offset+1):end] + 0.5 * rot * dsR * sh[(offset+1):end]
    @. b0[(offset+1):end] = @views b[(offset+1):end] + 0.5 * rot * dsR * sb[(offset+1):end]
    (h0, b0)
end
function calc_w0_2DF_2D(df0::AM, midpoint::AM, weight::AV)
    w = Vector{Float64}(undef, 4)
    u = @views midpoint[:, 1]
    v = @views midpoint[:, 2]
    h = @views df0[:, 1]
    b = @views df0[:, 2]
    w[1] = sum(@. weight * h)
    w[2] = sum(@. weight * u * h)
    w[3] = sum(@. weight * v * h)
    w[4] = 0.5 * sum(@. weight * ((u .^ 2 + v .^ 2) * h + b))
    w
end
function calc_w0_2DF_2D(h::AV, b::AV, u::AV, v::AV, weight::AV)
    w = Vector{Float64}(undef, 4)
    w[1] = sum(@. weight * h)
    w[2] = sum(@. weight * u * h)
    w[3] = sum(@. weight * v * h)
    w[4] = 0.5 * sum(@. weight * ((u .^ 2 + v .^ 2) * h + b))
    w
end
function calc_qf_2DF_2D(df0::AM, prim::AV, midpoint::AM, weight::AV)
    q = Vector{Float64}(undef, DIM)
    h = @views df0[:, 1]
    b = @views df0[:, 2]
    u = @views midpoint[:, 1]
    v = @views midpoint[:, 2]
    q[1] =
        0.5 * (
            sum(@. weight * (u - prim[2]) * ((u - prim[2])^2 + (v - prim[3])^2) * h) +
            sum(@. weight * (u - prim[2]) * b)
        )
    q[2] =
        0.5 * (
            sum(@. weight * (v - prim[3]) * ((u - prim[2])^2 + (v - prim[3])^2) * h) +
            sum(@. weight * (v - prim[3]) * b)
        )
    q
end
function calc_qf_2DF_2D(h::AV, b, prim::AV, u::AV, v, weight::AV)
    q = Vector{Float64}(undef, DIM)
    q[1] =
        0.5 * (
            sum(@. weight * (u - prim[2]) * ((u - prim[2])^2 + (v - prim[3])^2) * h) +
            sum(@. weight * (u - prim[2]) * b)
        )
    q[2] =
        0.5 * (
            sum(@. weight * (v - prim[3]) * ((u - prim[2])^2 + (v - prim[3])^2) * h) +
            sum(@. weight * (v - prim[3]) * b)
        )
    q
end
function calc_a(w0::AV, prim0::AV, cLw::AV, cRw, dxL::Real, dxR, K::Real, rot::Float64)
    sw = @. rot * (cLw - w0) / (0.5 * dxL)
    aL = micro_slope(sw, prim0, K)
    sw = @. rot * (w0 - cRw) / (0.5 * dxR)
    aR = micro_slope(sw, prim0, K)
    (aL, aR)
end
function micro_slope_2D(sw::AV, prim::AV, K::Real)
    a = Vector{Float64}(undef, 4)
    a[4] =
        4 * prim[4]^2 / (K + 2) / prim[1] * (
            2 * sw[4] - 2 * prim[2] * sw[2] - 2 * prim[3] * sw[3] +
            sw[1] * (prim[2]^2 + prim[3]^2 - 0.5 * (K + 2) / prim[4])
        )
    a[3] = 2 * prim[4] / prim[1] * (sw[3] - prim[3] * sw[1]) - prim[3] * a[4]
    a[2] = 2 * prim[4] / prim[1] * (sw[2] - prim[2] * sw[1]) - prim[2] * a[4]
    a[1] =
        sw[1] / prim[1] - prim[2] * a[2] - prim[3] * a[3] -
        0.5 * (prim[2]^2 + prim[3]^2 + 0.5 * (K + 2) / prim[4]) * a[4]
    a
end
function moment_u_2D(prim::AV, n::Int, m::Int, K::Float64, rot::Float64, dir::Int)# n>=1 calculate <u^n> <u^n>_>0 <u^n>_<0 <v^m> <ξ^4>
    if dir == 1
        Mu_L = Array{Float64}(undef, n + 1)
        Mu = Array{Float64}(undef, n + 1)
        Mu_R = Array{Float64}(undef, n + 1)
        Mu_L[1] = 0.5 * erfc(rot * √(prim[4]) * prim[2])
        Mu_L[2] =
            prim[2] * Mu_L[1] - rot * 0.5 * exp(-prim[4] * prim[2]^2) / (√(π * prim[4]))
        Mu_R[1] = 0.5 * erfc(-rot * √(prim[4]) * prim[2])
        Mu_R[2] =
            prim[2] * Mu_R[1] + rot * 0.5 * exp(-prim[4] * prim[2]^2) / (√(π * prim[4]))
        for i = 1:(n-1)
            Mu_L[i+2] = prim[2] * Mu_L[i+1] + 0.5 * i / prim[4] * Mu_L[i]
            Mu_R[i+2] = prim[2] * Mu_R[i+1] + 0.5 * i / prim[4] * Mu_R[i]
        end
        Mu = @. Mu_L + Mu_R
        Mv = Array{Float64}(undef, m + 1)
        Mv[1] = one(prim[3])
        Mv[2] = prim[3]
        for i = 1:(m-1)
            Mv[i+2] = prim[3] * Mv[i+1] + 0.5 * i * Mv[i] / prim[4]
        end
    else
        Mu = Array{Float64}(undef, m + 1)
        Mu_L = Array{Float64}(undef, n + 1)
        Mv = Array{Float64}(undef, n + 1)
        Mu_R = Array{Float64}(undef, n + 1)
        Mu_L[1] = 0.5 * erfc(rot * √(prim[4]) * prim[3])
        Mu_L[2] =
            prim[3] * Mu_L[1] - rot * 0.5 * exp(-prim[4] * prim[3]^2) / (√(π * prim[4]))
        Mu_R[1] = 0.5 * erfc(-rot * √(prim[4]) * prim[3])
        Mu_R[2] =
            prim[3] * Mu_R[1] + rot * 0.5 * exp(-prim[4] * prim[3]^2) / (√(π * prim[4]))
        for i = 1:(n-1)
            Mu_L[i+2] = prim[3] * Mu_L[i+1] + 0.5 * i / prim[4] * Mu_L[i]
            Mu_R[i+2] = prim[3] * Mu_R[i+1] + 0.5 * i / prim[4] * Mu_R[i]
        end
        Mv = @. Mu_L + Mu_R
        Mu[1] = one(prim[2])
        Mu[2] = prim[2]
        for i = 1:(m-1)
            Mu[i+2] = prim[2] * Mu[i+1] + 0.5 * i * Mu[i] / prim[4]
        end
    end

    Mξ2 = 0.5 * K / (prim[4])
    Mξ4 = (K^2 + 2.0 * K) / (4.0 * prim[4]^2)
    Mξ = [one(Mξ2), Mξ2, Mξ4]

    (Mu, Mv, Mξ, Mu_L, Mu_R)
end
function moment_uv(Mu, Mv, Mξ, α::Int, β::Int, δ::Int) # calc <u^αv^βξ^δψ>
    Muv = Vector{Float64}(undef, 4)
    Muv[1] = Mu[α+1] * Mv[β+1] * Mξ[Int(δ / 2 + 1)]
    Muv[2] = Mu[α+2] * Mv[β+1] * Mξ[Int(δ / 2 + 1)]
    Muv[3] = Mu[α+1] * Mv[β+2] * Mξ[Int(δ / 2 + 1)]
    Muv[4] =
        0.5 * (
            Mu[α+3] * Mv[β+1] * Mξ[Int(δ / 2 + 1)] +
            Mu[α+1] * Mv[β+3] * Mξ[Int(δ / 2 + 1)] +
            Mu[α+1] * Mv[β+1] * Mξ[Int(δ / 2 + 2)]
        )
    Muv
end
function moment_au(a, Mu, Mv, Mξ, α, β)# calc <au^αv^βψ> 
    a[1] * moment_uv(Mu, Mv, Mξ, α, β, 0) +
    a[2] * moment_uv(Mu, Mv, Mξ, α + 1, β, 0) +
    a[3] * moment_uv(Mu, Mv, Mξ, α, β + 1, 0) +
    0.5 * a[4] * moment_uv(Mu, Mv, Mξ, α + 2, β, 0) +
    0.5 * a[4] * moment_uv(Mu, Mv, Mξ, α, β + 2, 0) +
    0.5 * a[4] * moment_uv(Mu, Mv, Mξ, α, β, 2)
end
function calc_A(prim::AV, aL, aR, Mu, Mv, Mu_L, Mu_R, Mξ, K, dir::Int)
    if dir == 1
        Mau_L = moment_au(aL, Mu_L, Mv, Mξ, 1, 0)
        Mau_R = moment_au(aR, Mu_R, Mv, Mξ, 1, 0)
    else
        Mau_L = moment_au(aL, Mu, Mu_L, Mξ, 0, 1)
        Mau_R = moment_au(aR, Mu, Mu_R, Mξ, 0, 1)
    end
    sw = -prim[1] * (Mau_L + Mau_R)
    micro_slope(sw, prim, K)
end
function calc_time_int(τ, Δt)
    Mt = Vector{Float64}(undef, 5)
    Mt[4] = τ * (1.0 - exp(-Δt / τ))
    Mt[5] = -τ * Δt * exp(-Δt / τ) + τ * Mt[4]
    Mt[1] = Δt - Mt[4]
    Mt[2] = -τ * Mt[1] + Mt[5]
    Mt[3] = 0.5 * Δt^2 - τ * Mt[1]
    Mt
end
function calc_flux_g0(prim, Mt, Mu, Mu_L, Mu_R, Mv, Mξ, aL, aR, A, dir)
    if dir == 1
        Mau_0 = moment_uv(Mu, Mv, Mξ, 1, 0, 0)
        Mau2_L = moment_au(aL, Mu_L, Mv, Mξ, 2, 0)
        Mau2_R = moment_au(aR, Mu_R, Mv, Mξ, 2, 0)
        Mau_T = moment_au(A, Mu, Mv, Mξ, 1, 0)
    else
        Mau_0 = moment_uv(Mu, Mv, Mξ, 0, 1, 0)
        Mau2_L = moment_au(aL, Mu, Mu_L, Mξ, 0, 2)
        Mau2_R = moment_au(aR, Mu, Mu_R, Mξ, 0, 2)
        Mau_T = moment_au(A, Mu, Mv, Mξ, 0, 1)
    end
    Mt[1] * prim[1] * Mau_0 + Mt[2] * prim[1] * (Mau2_L + Mau2_R) + Mt[3] * prim[1] * Mau_T
end
function calc_flux_f0_2DF_2D(Mt, H⁺::AV, B⁺, u::AV, v, weight, h::AV, b, sh, sb::AV, vn::AV)
    F = Vector{Float64}(undef, 4)
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
        (sum(weight .* vn .^ 2 .* (u .^ 2 + v .^ 2) .* sh) + sum(weight .* vn .^ 2 .* sb))
    F
end

function calc_fhb(
    H0::AV,
    B0,
    H_plus,
    B_plus,
    Mt::AV,
    Mξ,
    u::AV,
    v,
    h::AV,
    b::AV,
    sh::AV,
    sb::AV,
    aL::AV,
    aR,
    aT::AV,
    offset::Int,
    vn::AV,
    ds::Float64,
)
    Θ = Array{Float64}(undef, length(H0))
    Θ[1:offset] .= 1.0
    Θ[(offset+1):end] .= 0.0
    Fh = @. Mt[1] * vn * (H0 + H_plus) +
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
    Fb = @. Mt[1] * vn * (B0 + B_plus) +
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
    (Fh .* ds, Fb .* ds)
end

function calc_time_int_O1(τ0::Float64, Δt::Float64)
    Mt = Vector{Float64}(undef, 2)
    Mt[2] = τ0 * (1.0 - exp(-Δt / τ0))
    Mt[1] = Δt - Mt[2]
    Mt
end
function calc_flux_g0_O1(prim, Mt, Mu, Mv, Mξ, dir)
    if dir == 1
        Mau_0 = moment_uv(Mu, Mv, Mξ, 1, 0, 0)
    else
        Mau_0 = moment_uv(Mu, Mv, Mξ, 0, 1, 0)
    end
    Mt[1] * prim[1] * Mau_0
end
function moment_u_O1(prim::AV, n::Int, m::Int, K::Float64, rot::Float64, dir::Int)# n>=1 calculate <u^n> <u^n>_>0 <u^n>_<0 <v^m> <ξ^4>
    if dir == 1
        Mu = Array{Float64}(undef, n + 1)
        Mu[1] = one(prim[2])
        Mu[2] = prim[2]
        for i = 1:(n-1)
            Mu[i+2] = prim[2] * Mu[i+1] + 0.5 * i * Mu[i] / prim[4]
        end
        Mv = Array{Float64}(undef, m + 1)
        Mv[1] = one(prim[3])
        Mv[2] = prim[3]
        for i = 1:(m-1)
            Mv[i+2] = prim[3] * Mv[i+1] + 0.5 * i * Mv[i] / prim[4]
        end
    else
        Mv = Array{Float64}(undef, n + 1)
        Mv[1] = one(prim[3])
        Mv[2] = prim[3]
        for i = 1:(n-1)
            Mv[i+2] = prim[3] * Mv[i+1] + 0.5 * i * Mv[i] / prim[4]
        end
        Mu = Array{Float64}(undef, m + 1)
        Mu[1] = one(prim[2])
        Mu[2] = prim[2]
        for i = 1:(m-1)
            Mu[i+2] = prim[2] * Mu[i+1] + 0.5 * i * Mu[i] / prim[4]
        end
    end

    Mξ2 = 0.5 * K / (prim[4])
    Mξ4 = (K^2 + 2.0 * K) / (4.0 * prim[4]^2)
    Mξ = [one(Mξ2), Mξ2, Mξ4]

    (Mu, Mv, Mξ)
end
function calc_flux_f0_O1(Mt, H⁺::AV, B⁺, u::AV, v, weight, h::AV, b, vn::AV)
    F = Vector{Float64}(undef, 4)
    F[1] = Mt[1] * sum(weight .* vn .* H⁺) + Mt[2] * sum(weight .* vn .* h)
    F[2] = Mt[1] * sum(weight .* vn .* u .* H⁺) + Mt[2] * sum(weight .* vn .* u .* h)
    F[3] = Mt[1] * sum(weight .* vn .* v .* H⁺) + Mt[2] * sum(weight .* vn .* v .* h)
    F[4] =
        Mt[1] *
        0.5 *
        (sum(weight .* vn .* (u .^ 2 + v .^ 2) .* H⁺) + sum(weight .* vn .* B⁺)) +
        Mt[2] * 0.5 * (sum(weight .* vn .* (u .^ 2 + v .^ 2) .* h) + sum(weight .* vn .* b))
    F
end
function calc_fhb_O1(
    H0::AV,
    B0,
    H_plus,
    B_plus,
    Mt::AV,
    h::AV,
    b::AV,
    offset::Int,
    vn::AV,
    ds::Float64,
)
    Θ = Array{Float64}(undef, length(H0))
    Θ[1:offset] .= 1.0
    Θ[(offset+1):end] .= 0.0
    Fh = @. Mt[1] * vn * (H0 + H_plus) + Mt[4] * vn * h
    Fb = @. Mt[1] * vn * (B0 + B_plus) + Mt[4] * vn * b
    (Fh .* ds, Fb .* ds)
end
