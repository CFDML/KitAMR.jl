using KitAMR
using LinearAlgebra
using Printf

function make_velocity_grid(n::Int; vmin=-5.0, vmax=5.0)
    du = (vmax - vmin) / n
    midpoint = Matrix{Float64}(undef, n * n, 2)
    weight = fill(du * du, n * n)
    row = 0
    for j in 1:n, i in 1:n
        row += 1
        midpoint[row, 1] = vmin + (i - 0.5) * du
        midpoint[row, 2] = vmin + (j - 0.5) * du
    end
    return midpoint, weight
end

function make_distribution(midpoint; prim=(1.0, 0.25, -0.15, 1.0), K=2.0)
    ρ, U, V, λ = prim
    df = Matrix{Float64}(undef, size(midpoint, 1), 2)
    @inbounds for i in axes(midpoint, 1)
        u = midpoint[i, 1]
        v = midpoint[i, 2]
        c1 = u - U
        c2 = v - V
        c_sq = c1 * c1 + c2 * c2
        h0 = ρ * λ / π * exp(-λ * c_sq)
        heat_skew = exp(0.035 * c1 * (c_sq - 2.0) - 0.02 * c2 * (c_sq - 2.0))
        df[i, 1] = h0 * heat_skew
        df[i, 2] = df[i, 1] * K / (2λ)
    end
    return df
end

function moments_2D2F(midpoint, df, weight)
    w = zeros(4)
    U = zeros(2)
    @inbounds for i in axes(df, 1)
        h = df[i, 1]
        b = df[i, 2]
        u = midpoint[i, 1]
        v = midpoint[i, 2]
        ω = weight[i]
        w[1] += ω * h
        w[2] += ω * u * h
        w[3] += ω * v * h
        w[4] += 0.5 * ω * ((u * u + v * v) * h + b)
    end
    U[1] = w[2] / w[1]
    U[2] = w[3] / w[1]
    qf = zeros(2)
    @inbounds for i in axes(df, 1)
        h = df[i, 1]
        b = df[i, 2]
        c1 = midpoint[i, 1] - U[1]
        c2 = midpoint[i, 2] - U[2]
        c_sq = c1 * c1 + c2 * c2
        ω = weight[i]
        qf[1] += 0.5 * ω * c1 * (c_sq * h + b)
        qf[2] += 0.5 * ω * c2 * (c_sq * h + b)
    end
    return w, qf
end

function make_vsdata(midpoint, weight, df)
    n = size(df, 1)
    return VsData{2,2}(
        n,
        zeros(Int8, n),
        copy(weight),
        copy(midpoint),
        copy(df),
        zeros(n, 2, 2),
        zeros(n, 2),
    )
end

function bench!(vs_data, base_df, w, qf; heat_constraint::Bool, reps::Int)
    t0 = time_ns()
    for _ in 1:reps
        vs_data.df .= base_df
        if heat_constraint
            KitAMR.conserved_I_porjection!(vs_data, w, qf)
        else
            KitAMR.conserved_I_porjection!(vs_data, w)
        end
    end
    return (time_ns() - t0) / reps / 1e6
end

function main()
    n = parse(Int, get(ENV, "KITAMR_IPROJ_BENCH_N", "24"))
    reps = parse(Int, get(ENV, "KITAMR_IPROJ_BENCH_REPS", "40"))
    midpoint, weight = make_velocity_grid(n)
    target_df = make_distribution(midpoint)
    target_w, target_qf = moments_2D2F(midpoint, target_df, weight)

    mapped_df = copy(target_df)
    @inbounds for i in axes(mapped_df, 1)
        u = midpoint[i, 1]
        v = midpoint[i, 2]
        mapped_df[i, 1] *= exp(0.03 * sin(1.7u) * cos(1.3v))
    end

    vs_old = make_vsdata(midpoint, weight, mapped_df)
    vs_new = make_vsdata(midpoint, weight, mapped_df)

    bench!(vs_old, mapped_df, target_w, target_qf; heat_constraint=false, reps=3)
    bench!(vs_new, mapped_df, target_w, target_qf; heat_constraint=true, reps=3)

    t_old = bench!(vs_old, mapped_df, target_w, target_qf; heat_constraint=false, reps=reps)
    t_new = bench!(vs_new, mapped_df, target_w, target_qf; heat_constraint=true, reps=reps)

    old_w, old_qf = moments_2D2F(midpoint, vs_old.df, weight)
    new_w, new_qf = moments_2D2F(midpoint, vs_new.df, weight)

    @printf("Velocity grid: %d x %d (%d cells), reps=%d\n", n, n, n * n, reps)
    @printf("No heat constraint:   %.4f ms/solve, |Δw|=%.3e, |Δq|=%.3e\n",
        t_old, norm(old_w - target_w), norm(old_qf - target_qf))
    @printf("With heat constraint: %.4f ms/solve, |Δw|=%.3e, |Δq|=%.3e\n",
        t_new, norm(new_w - target_w), norm(new_qf - target_qf))
    @printf("Cost ratio: %.3fx\n", t_new / t_old)
end

main()
