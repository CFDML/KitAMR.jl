"""
    build_psi_matrix(midpoint)

构造 Ψ ∈ ℝ^{N×(D+2)}, 第 k 行为 ψ_k = (1, ξ_k, |ξ_k|²/2).
"""
function build_psi_matrix(midpoint::AbstractMatrix{T}) where {T<:Real}
    N, D = size(midpoint)
    Psi = Matrix{T}(undef, N, D + 2)
    @inbounds for k = 1:N
        Psi[k, 1] = one(T)
        s = zero(T)
        for d = 1:D
            x = midpoint[k, d]
            Psi[k, 1+d] = x
            s += x * x
        end
        Psi[k, D+2] = s / 2
    end
    return Psi
end

"""
    dual_objective(Psi, f, W, weight, lambda)

I-投影的对偶目标
Φ(λ) = ⟨f* exp(λ·ψ)⟩ - λ·W = Σ_k ω_k f_k exp(λ·ψ_k) - λ·W.

满足 ∇Φ = G, ∇²Φ = J. 严格凸, Newton 方向是 Φ 的下降方向, 可用于 Armijo.
"""
function dual_objective(
    Psi::AbstractMatrix{T},
    f::AbstractVector{T},
    W::AbstractVector{T},
    weight::AbstractVector{T},
    lambda::AbstractVector{T},
) where {T<:Real}
    N, M = size(Psi)
    acc = zero(T)
    @inbounds for k = 1:N
        dot_lp = zero(T)
        for i = 1:M
            dot_lp += lambda[i] * Psi[k, i]
        end
        acc += weight[k] * f[k] * exp(dot_lp)
    end
    return acc - dot(lambda, W)
end

"""
    solve_iprojection(midpoint, f, W, weight; tol=1e-10, maxiter=20)

Newton 迭代求解 I-投影对偶方程 G(λ) = ⟨ψ f exp(λ·ψ)⟩ - W = 0, 初值 λ=0.
返回收敛后的 λ ∈ ℝ^{D+2}.
"""
function solve_I_projection(
    midpoint::AbstractMatrix{T},
    f::AbstractVector{T},
    W::AbstractVector{T},
    weight::AbstractVector{T};
    tol::Real = 1e-10,
    maxiter::Integer = 10,
) where {T<:Real}
    N, D = size(midpoint)
    M = D + 2
    Psi = build_psi_matrix(midpoint)
    lambda = zeros(T, M)
    G = Vector{T}(undef, M)
    J = Matrix{T}(undef, M, M)
    tol_eff = tol * max(one(T), norm(W))
    G_prev_norm = T(Inf)
    stall_count = 0

    # shave negative distribution
    f_min = 1.1*minimum(f)
    if f_min<0
        fp_min = minimum(x for x in f if x>0)
        d = fp_min-f_min
        for i in eachindex(f)
            if f[i]<0
                f[i] = (f[i]-f_min)/d*fp_min
            end
        end
    end


    for _ = 1:maxiter
        fill!(G, zero(T))
        fill!(J, zero(T))

        @inbounds for k = 1:N
            dot_lp = zero(T)
            for i = 1:M
                dot_lp += lambda[i] * Psi[k, i]
            end
            c = weight[k] * f[k] * exp(dot_lp)

            for i = 1:M
                psi_ki = Psi[k, i]
                ci = c * psi_ki
                G[i] += ci
                for j = i:M
                    J[i, j] += ci * Psi[k, j]
                end
            end
        end

        G .-= W
        G_norm = norm(G)

        if G_norm < tol_eff
            return lambda, Psi
        end

        # 检测停滞: 残差不再显著下降
        if G_norm > 0.9 * G_prev_norm
            stall_count += 1
            if stall_count ≥ 2
                return lambda, Psi    # 已达精度极限, 提前退出
            end
        else
            stall_count = 0
        end
        G_prev_norm = G_norm
        Δλ = -(Symmetric(J, :U) \ G)
        # Armijo 回溯线搜索
        Φ0 = dual_objective(Psi, f, W, weight, lambda)
        slope = dot(G, Δλ)              # = ∇Φ·Δλ < 0
        α = one(T)
        for _ = 1:maxiter
            if dual_objective(Psi, f, W, weight, lambda + α*Δλ) ≤ Φ0 + 1e-4*α*slope
                break
            end
            α *= 0.5
        end
        lambda .+= α * Δλ
    end
    if norm(G) > 1e-6
        @warn "Newton 未在 $maxiter 步内收敛, 残差 $(norm(G))"
    end
    return lambda, Psi
end


function conserved_I_porjection!(vs_data::VsData{2,2}, w::AbstractVector)
    e_int = @views dot(vs_data.weight, vs_data.df[:, 2])/2
    w_trans = copy(w);
    w_trans[end]-=e_int
    λ, Ψ = @views solve_I_projection(
        vs_data.midpoint,
        vs_data.df[:, 1],
        w_trans,
        vs_data.weight,
    )
    for i in axes(vs_data.df, 1)
        vs_data.df[i, 1] *= @views exp(dot(λ, Ψ[i, :]))
    end
end

function conserved_I_porjection!(vs_data::VsData{3,1}, w::AbstractVector)
    λ, Ψ = @views solve_I_projection(vs_data.midpoint, vs_data.df[:, 1], w, vs_data.weight)
    for i in axes(vs_data.df, 1)
        vs_data.df[i, 1] *= @views exp(dot(λ, Ψ[i, :]))
    end
end


function iterate!(::Type{CIP_Marching}, ka::KA)
    kinfo = ka.kinfo
    gas = kinfo.config.gas
    trees = ka.kdata.field.trees
    Δt = kinfo.status.Δt
    for i in eachindex(trees.data)
        @inbounds for j in eachindex(trees.data[i])
            ps_data = trees.data[i][j]
            isa(ps_data, InsideSolidData) && continue
            ps_data.bound_enc<0 && continue
            vs_data = ps_data.vs_data
            area = reduce(*, ps_data.ds)
            ps_data.w .+= ps_data.flux .* Δt / area # Macroscopic update
            positivity_preserving_ib!(ps_data, area, Δt)
            prim_c = get_prim(ps_data, kinfo) # Conserved macroscopic variables
            f = vs_data.df
            f .+= Δt/area*vs_data.flux # Convection first
            F_c = discrete_maxwell(vs_data.midpoint, prim_c, kinfo)
            ps_data.qf .= qf = calc_qf(vs_data, prim_c) # Heatflux after convection
            F_c .+= shakhov_part(vs_data.midpoint, F_c, prim_c, qf, kinfo) # g^{S,n+1}
            conserved_I_porjection!(vs_data, ps_data.w)
            # Collision process
            τ = get_τ(prim_c, gas.μᵣ, gas.ω) # τ^{n+1}
            f .*= τ/(τ+Δt)
            @. f += Δt/(τ+Δt)*F_c
            residual_check!(ps_data, prim_c, kinfo)
            ps_data.prim .= prim_c
            ps_data.flux .= 0.0
            vs_data.flux .= 0.0
        end
    end
    return nothing
end
