# Post-process 2D/3D vacuum-jet AMR efficiency runs.
using JLD2
using Printf

const OUTDIR = isempty(ARGS) ? joinpath(@__DIR__, "out") : ARGS[1]
const DEFAULT_MODE_ORDER = ["both_on", "both_off", "both_on_3d", "both_off_3d"]

metric(m, name::Symbol, default) = hasproperty(m, name) ? getproperty(m, name) : default
run_dim(m) = Int(metric(m, :dimension, length(metric(m, :geometry, Float64[])) ÷ 2))
field_cols(dim::Integer) = 3 * dim + 2
rho_col(dim::Integer) = 2 * dim + 1
theta_col(dim::Integer) = 3 * dim + 2

function discover_labels()
    explicit = strip(get(ENV, "VJE_ANALYZE_MODES", ""))
    if !isempty(explicit)
        return String.(split(replace(explicit, "," => " ")))
    end

    labels = String[]
    for label in DEFAULT_MODE_ORDER
        isfile(joinpath(OUTDIR, label, "metrics.jld2")) && push!(labels, label)
    end
    if isdir(OUTDIR)
        extras = sort([f for f in readdir(OUTDIR)
                       if isfile(joinpath(OUTDIR, f, "metrics.jld2"))])
        for label in extras
            label in labels || push!(labels, label)
        end
    end
    return labels
end

function load_mode(label)
    dir = joinpath(OUTDIR, label)
    mfile = joinpath(dir, "metrics.jld2")
    isfile(mfile) || return nothing
    m = load(mfile, "metrics")
    dim = run_dim(m)
    parts = sort([joinpath(dir, f) for f in readdir(dir)
                  if startswith(f, "field_rank") && endswith(f, ".jld2")])
    rows = isempty(parts) ? zeros(0, field_cols(dim)) :
        reduce(vcat, (load(p, "rows") for p in parts))
    if size(rows, 2) != field_cols(dim)
        @warn "unexpected field row width" label dim width=size(rows, 2) expected=field_cols(dim)
    end
    return (; label, m, rows, dim)
end

labels = discover_labels()
runs = Dict{String,Any}()
for label in labels
    r = load_mode(label)
    r === nothing || (runs[label] = r)
end
isempty(runs) && error("no runs found under $OUTDIR")

function reference_grid(run)
    m = run.m
    dim = run.dim
    geom = m.geometry
    nx = Int(m.ref_ps_x)
    ny = Int(m.ref_ps_y)
    x0, x1, y0, y1 = geom[1], geom[2], geom[3], geom[4]
    dx = (x1 - x0) / nx
    dy = (y1 - y0) / ny
    rcol = rho_col(dim)
    tcol = theta_col(dim)

    if dim == 2
        rho = fill(NaN, nx, ny)
        theta = fill(NaN, nx, ny)
        for k in axes(run.rows, 1)
            x = run.rows[k, 1]
            y = run.rows[k, 2]
            i = clamp(floor(Int, (x - x0) / dx) + 1, 1, nx)
            j = clamp(floor(Int, (y - y0) / dy) + 1, 1, ny)
            rho[i, j] = run.rows[k, rcol]
            theta[i, j] = run.rows[k, tcol]
        end
        return (; dim, nx, ny, x0, y0, dx, dy, rho, theta)
    elseif dim == 3
        nz = Int(metric(m, :ref_ps_z, 1))
        z0, z1 = geom[5], geom[6]
        dz = (z1 - z0) / nz
        rho = fill(NaN, nx, ny, nz)
        theta = fill(NaN, nx, ny, nz)
        for k in axes(run.rows, 1)
            x = run.rows[k, 1]
            y = run.rows[k, 2]
            z = run.rows[k, 3]
            i = clamp(floor(Int, (x - x0) / dx) + 1, 1, nx)
            j = clamp(floor(Int, (y - y0) / dy) + 1, 1, ny)
            l = clamp(floor(Int, (z - z0) / dz) + 1, 1, nz)
            rho[i, j, l] = run.rows[k, rcol]
            theta[i, j, l] = run.rows[k, tcol]
        end
        return (; dim, nx, ny, nz, x0, y0, z0, dx, dy, dz, rho, theta)
    end
    error("unsupported dimension $dim")
end

function is_reference_run(run)
    mode = String(metric(run.m, :mode, run.label))
    return replace(mode, "_3d" => "") == "both_off"
end

ref_by_dim = Dict{Int,Any}()
base_by_dim = Dict{Int,Any}()
for label in labels
    haskey(runs, label) || continue
    run = runs[label]
    if is_reference_run(run) && !haskey(base_by_dim, run.dim)
        base_by_dim[run.dim] = run.m
        ref_by_dim[run.dim] = reference_grid(run)
    end
end
for dim in sort(unique([run.dim for run in values(runs)]))
    if !haskey(ref_by_dim, dim)
        @warn "uniform reference is missing; accuracy and saving columns will be NaN" dim
    end
end

function ref_block_avg(g, x, y, hx, hy)
    i_lo = clamp(floor(Int, (x - hx - g.x0) / g.dx) + 1, 1, g.nx)
    i_hi = clamp(ceil(Int, (x + hx - g.x0) / g.dx), 1, g.nx)
    j_lo = clamp(floor(Int, (y - hy - g.y0) / g.dy) + 1, 1, g.ny)
    j_hi = clamp(ceil(Int, (y + hy - g.y0) / g.dy), 1, g.ny)
    sr = 0.0
    st = 0.0
    n = 0
    @inbounds for i in i_lo:i_hi, j in j_lo:j_hi
        rr = g.rho[i, j]
        tt = g.theta[i, j]
        (isnan(rr) || isnan(tt)) && continue
        sr += rr
        st += tt
        n += 1
    end
    n == 0 && return (NaN, NaN)
    return (sr / n, st / n)
end

function ref_block_avg(g, x, y, z, hx, hy, hz)
    i_lo = clamp(floor(Int, (x - hx - g.x0) / g.dx) + 1, 1, g.nx)
    i_hi = clamp(ceil(Int, (x + hx - g.x0) / g.dx), 1, g.nx)
    j_lo = clamp(floor(Int, (y - hy - g.y0) / g.dy) + 1, 1, g.ny)
    j_hi = clamp(ceil(Int, (y + hy - g.y0) / g.dy), 1, g.ny)
    k_lo = clamp(floor(Int, (z - hz - g.z0) / g.dz) + 1, 1, g.nz)
    k_hi = clamp(ceil(Int, (z + hz - g.z0) / g.dz), 1, g.nz)
    sr = 0.0
    st = 0.0
    n = 0
    @inbounds for i in i_lo:i_hi, j in j_lo:j_hi, k in k_lo:k_hi
        rr = g.rho[i, j, k]
        tt = g.theta[i, j, k]
        (isnan(rr) || isnan(tt)) && continue
        sr += rr
        st += tt
        n += 1
    end
    n == 0 && return (NaN, NaN)
    return (sr / n, st / n)
end

function errors_vs_ref(rows, g, dim)
    (g === nothing || g.dim != dim) && return (NaN, NaN, NaN, NaN)
    measure_sum = 0.0
    l1_rho = 0.0
    l2_rho = 0.0
    l1_theta = 0.0
    l2_theta = 0.0
    norm_rho = 0.0
    norm_theta = 0.0
    rcol = rho_col(dim)
    tcol = theta_col(dim)

    for k in axes(rows, 1)
        if dim == 2
            x, y = rows[k, 1], rows[k, 2]
            dx, dy = rows[k, 3], rows[k, 4]
            measure = dx * dy
            rr, tt = ref_block_avg(g, x, y, 0.5 * dx, 0.5 * dy)
        else
            x, y, z = rows[k, 1], rows[k, 2], rows[k, 3]
            dx, dy, dz = rows[k, 4], rows[k, 5], rows[k, 6]
            measure = dx * dy * dz
            rr, tt = ref_block_avg(g, x, y, z, 0.5 * dx, 0.5 * dy, 0.5 * dz)
        end
        (isnan(rr) || isnan(tt)) && continue
        dr = rows[k, rcol] - rr
        dt = rows[k, tcol] - tt
        measure_sum += measure
        l1_rho += abs(dr) * measure
        l2_rho += dr^2 * measure
        l1_theta += abs(dt) * measure
        l2_theta += dt^2 * measure
        norm_rho += abs(rr) * measure
        norm_theta += abs(tt) * measure
    end
    (measure_sum == 0 || norm_rho == 0 || norm_theta == 0) &&
        return (NaN, NaN, NaN, NaN)
    return (
        l1_rho / norm_rho,
        sqrt(l2_rho / measure_sum) / (norm_rho / measure_sum),
        l1_theta / norm_theta,
        sqrt(l2_theta / measure_sum) / (norm_theta / measure_sum),
    )
end

records = []
for label in labels
    haskey(runs, label) || continue
    run = runs[label]
    ref = get(ref_by_dim, run.dim, nothing)
    l1r, l2r, l1t, l2t = errors_vs_ref(run.rows, ref, run.dim)
    push!(records, (; label, run.dim, m = run.m, l1r, l2r, l1t, l2t))
end

println("\n", "="^118)
println(" VACUUM-JET AMR EFFICIENCY REPORT  OUTDIR=$OUTDIR")
println("="^118)
@printf("%-13s %3s %5s %6s %10s %10s %12s %10s %9s %9s %9s %9s\n",
        "mode", "dim", "np", "steps", "wall_s", "core_h", "phase", "PS",
        "dataGB", "dataMax", "memGB", "memMax")
println("-"^118)
for r in records
    m = r.m
    @printf("%-13s %3d %5d %6d %10.3f %10.5f %12d %10d %9.3f %9.3f %9.3f %9.3f\n",
            String(metric(m, :mode, r.label)), r.dim, m.np, m.n_steps, m.wall_s, m.core_h,
            m.total_phase, m.total_ps, m.data_sum_GB, m.data_max_GB, m.mem_sum_GB,
            m.mem_max_GB)
end
println("-"^118)
@printf("%-13s %3s %10s %10s %10s %10s %12s %12s %12s\n",
        "mode", "dim", "L1(rho)", "L2(rho)", "L1(theta)", "L2(theta)",
        "phase_save", "coreh_save", "data_save")
println("-"^118)
for r in records
    m = r.m
    base = get(base_by_dim, r.dim, nothing)
    phase_save = base === nothing ? NaN : base.total_phase / m.total_phase
    coreh_save = base === nothing ? NaN : base.core_h / max(m.core_h, eps())
    data_save = base === nothing ? NaN : base.data_sum_GB / max(m.data_sum_GB, eps())
    @printf("%-13s %3d %10.3e %10.3e %10.3e %10.3e %11.2fx %11.2fx %11.2fx\n",
            String(metric(m, :mode, r.label)), r.dim, r.l1r, r.l2r, r.l1t, r.l2t,
            phase_save, coreh_save, data_save)
end
println("="^118)

mkpath(OUTDIR)
csv = joinpath(OUTDIR, "amr_effi_report.csv")
open(csv, "w") do io
    println(io, "mode,label,dimension,np,steps,sim_time,wall_s,core_h,total_phase,total_ps,max_vs,data_sum_GB,data_max_GB,mem_sum_GB,mem_max_GB,L1_rho,L2_rho,L1_theta,L2_theta,phase_save,coreh_save,data_save,mass,energy,minrho,maxrho,maxspeed")
    for r in records
        m = r.m
        base = get(base_by_dim, r.dim, nothing)
        phase_save = base === nothing ? NaN : base.total_phase / m.total_phase
        coreh_save = base === nothing ? NaN : base.core_h / max(m.core_h, eps())
        data_save = base === nothing ? NaN : base.data_sum_GB / max(m.data_sum_GB, eps())
        @printf(io, "%s,%s,%d,%d,%d,%.8g,%.6f,%.8f,%d,%d,%d,%.8f,%.8f,%.8f,%.8f,%.8e,%.8e,%.8e,%.8e,%.6f,%.6f,%.6f,%.8e,%.8e,%.8e,%.8e,%.8e\n",
                String(metric(m, :mode, r.label)), String(metric(m, :label, r.label)), r.dim,
                m.np, m.n_steps, m.sim_time, m.wall_s, m.core_h, m.total_phase,
                m.total_ps, m.max_vs, m.data_sum_GB, m.data_max_GB, m.mem_sum_GB,
                m.mem_max_GB, r.l1r, r.l2r, r.l1t, r.l2t, phase_save, coreh_save,
                data_save, m.mass, m.energy, m.minrho, m.maxrho, m.maxspeed)
    end
end

tsv = joinpath(OUTDIR, "summary.tsv")
open(tsv, "w") do io
    @printf(io, "mode\tdimension\tnp\tsteps\twall_s\tcore_h\tPS\tphase\tmaxVS\tdata_sum_GB\tdata_max_GB\tmem_sum_GB\tmem_max_GB\n")
    for r in records
        m = r.m
        @printf(io, "%s\t%d\t%d\t%d\t%.6f\t%.8f\t%d\t%d\t%d\t%.8f\t%.8f\t%.8f\t%.8f\n",
                String(metric(m, :mode, r.label)), r.dim, m.np, m.n_steps, m.wall_s,
                m.core_h, m.total_ps, m.total_phase, m.max_vs, m.data_sum_GB,
                m.data_max_GB, m.mem_sum_GB, m.mem_max_GB)
    end
end

println("wrote $csv")
println("wrote $tsv")
