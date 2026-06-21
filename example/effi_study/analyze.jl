# Post-process R2D vacuum-jet AMR efficiency runs.
using JLD2
using Printf

const OUTDIR = isempty(ARGS) ? joinpath(@__DIR__, "out") : ARGS[1]
const MODES = ["both_on", "ps_only", "vs_only", "both_off"]

function load_mode(dir)
    mfile = joinpath(dir, "metrics.jld2")
    isfile(mfile) || return nothing
    m = load(mfile, "metrics")
    parts = sort([joinpath(dir, f) for f in readdir(dir)
                  if startswith(f, "field_rank") && endswith(f, ".jld2")])
    rows = isempty(parts) ? zeros(0, 8) : reduce(vcat, (load(p, "rows") for p in parts))
    return (; m, rows)
end

runs = Dict{String,Any}()
for mode in MODES
    r = load_mode(joinpath(OUTDIR, mode))
    r === nothing || (runs[mode] = r)
end
isempty(runs) && error("no runs found under $OUTDIR")

function reference_grid(run)
    m = run.m
    geom = m.geometry
    nx = m.ref_ps_x
    ny = m.ref_ps_y
    x0, x1, y0, y1 = geom[1], geom[2], geom[3], geom[4]
    dx = (x1 - x0) / nx
    dy = (y1 - y0) / ny
    rho = fill(NaN, nx, ny)
    theta = fill(NaN, nx, ny)
    for k in axes(run.rows, 1)
        x = run.rows[k, 1]
        y = run.rows[k, 2]
        i = clamp(floor(Int, (x - x0) / dx) + 1, 1, nx)
        j = clamp(floor(Int, (y - y0) / dy) + 1, 1, ny)
        rho[i, j] = run.rows[k, 5]
        theta[i, j] = run.rows[k, 8]
    end
    return (; nx, ny, x0, y0, dx, dy, rho, theta)
end

ref = haskey(runs, "both_off") ? reference_grid(runs["both_off"]) : nothing
ref === nothing && @warn "both_off reference is missing; accuracy columns will be NaN"

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

function errors_vs_ref(rows, g)
    g === nothing && return (NaN, NaN, NaN, NaN)
    area_sum = 0.0
    l1_rho = 0.0
    l2_rho = 0.0
    l1_theta = 0.0
    l2_theta = 0.0
    norm_rho = 0.0
    norm_theta = 0.0
    for k in axes(rows, 1)
        x = rows[k, 1]
        y = rows[k, 2]
        dx = rows[k, 3]
        dy = rows[k, 4]
        area = dx * dy
        rr, tt = ref_block_avg(g, x, y, 0.5 * dx, 0.5 * dy)
        (isnan(rr) || isnan(tt)) && continue
        dr = rows[k, 5] - rr
        dt = rows[k, 8] - tt
        area_sum += area
        l1_rho += abs(dr) * area
        l2_rho += dr^2 * area
        l1_theta += abs(dt) * area
        l2_theta += dt^2 * area
        norm_rho += abs(rr) * area
        norm_theta += abs(tt) * area
    end
    area_sum == 0 && return (NaN, NaN, NaN, NaN)
    return (
        l1_rho / norm_rho,
        sqrt(l2_rho / area_sum) / (norm_rho / area_sum),
        l1_theta / norm_theta,
        sqrt(l2_theta / area_sum) / (norm_theta / area_sum),
    )
end

records = []
for mode in MODES
    haskey(runs, mode) || continue
    m = runs[mode].m
    l1r, l2r, l1t, l2t = errors_vs_ref(runs[mode].rows, ref)
    push!(records, (; mode, m, l1r, l2r, l1t, l2t))
end

base = haskey(runs, "both_off") ? runs["both_off"].m : nothing

println("\n", "="^110)
println(" R2D VACUUM-JET AMR EFFICIENCY REPORT  OUTDIR=$OUTDIR")
println("="^110)
@printf("%-9s %5s %6s %10s %10s %12s %10s %9s %9s %9s %9s\n",
        "mode", "np", "steps", "wall_s", "core_h", "phase", "PS",
        "dataGB", "dataMax", "memGB", "memMax")
println("-"^110)
for r in records
    m = r.m
    @printf("%-9s %5d %6d %10.3f %10.5f %12d %10d %9.3f %9.3f %9.3f %9.3f\n",
            m.mode, m.np, m.n_steps, m.wall_s, m.core_h, m.total_phase, m.total_ps,
            m.data_sum_GB, m.data_max_GB, m.mem_sum_GB, m.mem_max_GB)
end
println("-"^110)
@printf("%-9s %10s %10s %10s %10s %12s %12s %12s\n",
        "mode", "L1(rho)", "L2(rho)", "L1(theta)", "L2(theta)",
        "phase_save", "coreh_save", "data_save")
println("-"^110)
for r in records
    m = r.m
    phase_save = base === nothing ? NaN : base.total_phase / m.total_phase
    coreh_save = base === nothing ? NaN : base.core_h / max(m.core_h, eps())
    data_save = base === nothing ? NaN : base.data_sum_GB / max(m.data_sum_GB, eps())
    @printf("%-9s %10.3e %10.3e %10.3e %10.3e %11.2fx %11.2fx %11.2fx\n",
            m.mode, r.l1r, r.l2r, r.l1t, r.l2t, phase_save, coreh_save, data_save)
end
println("="^110)

mkpath(OUTDIR)
csv = joinpath(OUTDIR, "amr_effi_report.csv")
open(csv, "w") do io
    println(io, "mode,label,np,steps,sim_time,wall_s,core_h,total_phase,total_ps,max_vs,data_sum_GB,data_max_GB,mem_sum_GB,mem_max_GB,L1_rho,L2_rho,L1_theta,L2_theta,phase_save,coreh_save,data_save,mass,energy,minrho,maxrho,maxspeed")
    for r in records
        m = r.m
        phase_save = base === nothing ? NaN : base.total_phase / m.total_phase
        coreh_save = base === nothing ? NaN : base.core_h / max(m.core_h, eps())
        data_save = base === nothing ? NaN : base.data_sum_GB / max(m.data_sum_GB, eps())
        @printf(io, "%s,%s,%d,%d,%.8g,%.6f,%.8f,%d,%d,%d,%.8f,%.8f,%.8f,%.8f,%.8e,%.8e,%.8e,%.8e,%.6f,%.6f,%.6f,%.8e,%.8e,%.8e,%.8e,%.8e\n",
                m.mode, m.label, m.np, m.n_steps, m.sim_time, m.wall_s, m.core_h,
                m.total_phase, m.total_ps, m.max_vs,
                m.data_sum_GB, m.data_max_GB, m.mem_sum_GB, m.mem_max_GB,
                r.l1r, r.l2r, r.l1t, r.l2t, phase_save, coreh_save, data_save,
                m.mass, m.energy, m.minrho, m.maxrho, m.maxspeed)
    end
end

tsv = joinpath(OUTDIR, "summary.tsv")
open(tsv, "w") do io
    @printf(io, "mode\tnp\tsteps\twall_s\tcore_h\tPS\tphase\tmaxVS\tdata_sum_GB\tdata_max_GB\tmem_sum_GB\tmem_max_GB\n")
    for r in records
        m = r.m
        @printf(io, "%s\t%d\t%d\t%.6f\t%.8f\t%d\t%d\t%d\t%.8f\t%.8f\t%.8f\t%.8f\n",
                m.mode, m.np, m.n_steps, m.wall_s, m.core_h, m.total_ps,
                m.total_phase, m.max_vs, m.data_sum_GB, m.data_max_GB,
                m.mem_sum_GB, m.mem_max_GB)
    end
end

println("wrote $csv")
println("wrote $tsv")
