# Automatic AMR intervals use a short cached startup value from `Status`; after an AMR pass
# this file refreshes the cache from kinetic transport statistics and may relax it.
const AUTO_AMR_MIN_INTERVAL = 1
const AUTO_AMR_RATE_BINS_PER_OCTAVE = 8
const AUTO_AMR_RATE_HIST_OCTAVES = 32
const AUTO_AMR_RATE_HIST_BINS = AUTO_AMR_RATE_BINS_PER_OCTAVE * AUTO_AMR_RATE_HIST_OCTAVES

function _positive_interval(value::Integer, name::Symbol)
    ivalue = Int(value)
    ivalue > 0 || error("`$name` must be positive; got $ivalue.")
    return ivalue
end
function _positive_interval(value, name::Symbol)
    error("`$name` must be a positive integer, `:auto`, or a function returning a positive integer; got $(repr(value)).")
end

function _call_interval_function(fn::Function, p4est::P_pxest_t, ka::KA, kind::Symbol)
    if applicable(fn, p4est, ka, kind)
        return fn(p4est, ka, kind)
    elseif applicable(fn, p4est, ka)
        return fn(p4est, ka)
    elseif applicable(fn, ka, kind)
        return fn(ka, kind)
    elseif applicable(fn, ka)
        return fn(ka)
    end
    error("AMR interval function must accept one of `(p4est, ka, kind)`, `(p4est, ka)`, `(ka, kind)`, or `(ka)`.")
end

@inline _interval_name(kind::Symbol) = kind === :ps ? :ps_interval : :vs_interval
@inline _cached_auto_interval(ka::KA, kind::Symbol) =
    kind === :ps ? ka.kinfo.status.ps_interval_cached : ka.kinfo.status.vs_interval_cached

function _interval_value(interval::Integer, p4est::P_pxest_t, ka::KA, kind::Symbol)
    return _positive_interval(interval, _interval_name(kind))
end
function _interval_value(interval::Symbol, p4est::P_pxest_t, ka::KA, kind::Symbol)
    interval === :auto && return _cached_auto_interval(ka, kind)
    error("`$(kind)_interval` must be a positive integer, `:auto`, or a function; got $(repr(interval)).")
end
function _interval_value(fn::Function, p4est::P_pxest_t, ka::KA, kind::Symbol)
    value = _call_interval_function(fn, p4est, ka, kind)
    return _interval_function_value(value, p4est, ka, kind)
end
function _interval_value(interval, p4est::P_pxest_t, ka::KA, kind::Symbol)
    error("`$(kind)_interval` must be a positive integer, `:auto`, or a function; got $(repr(interval)).")
end
function _interval_function_value(value::Integer, p4est::P_pxest_t, ka::KA, kind::Symbol)
    return _positive_interval(value, _interval_name(kind))
end
function _interval_function_value(value::Symbol, p4est::P_pxest_t, ka::KA, kind::Symbol)
    return _interval_value(value, p4est, ka, kind)
end
function _interval_function_value(value, p4est::P_pxest_t, ka::KA, kind::Symbol)
    error("`$(kind)_interval` function must return a positive integer or `:auto`; got $(repr(value)).")
end

function _partition_interval_value(interval::Integer, ps_interval::Integer)
    return _positive_interval(interval, :partition_interval)
end
function _partition_interval_value(interval::Symbol, ps_interval::Integer)
    interval === :auto && return 2 * ps_interval
    error("`partition_interval` must be a positive integer or `:auto`; got $(repr(interval)).")
end
function _partition_interval_value(interval, ps_interval::Integer)
    error("`partition_interval` must be a positive integer or `:auto`; got $(repr(interval)).")
end

@inline _auto_amr_uses_default_pair(ps_interval, vs_interval) =
    ps_interval === :auto && vs_interval === :auto

@inline function _auto_amr_should_align_intervals(ka::KA, ps_interval, vs_interval)
    solver = ka.kinfo.config.solver
    return solver.AUTO_AMR_ALIGN_INTERVALS &&
           solver.PS_DYNAMIC_AMR &&
           solver.VS_DYNAMIC_AMR &&
           _auto_amr_uses_default_pair(ps_interval, vs_interval)
end

@inline function _auto_amr_multiple_pair(ps_interval::Integer, vs_interval::Integer)
    # Keep the two automatic estimates distinct, but make the longer one an integer multiple of
    # the shorter one.  Use the lower multiple so the transport-based interval estimate is never
    # relaxed by alignment; the scheduler can then co-trigger the two AMR stages without adding
    # standalone recovery passes.
    ps = Int(ps_interval)
    vs = Int(vs_interval)
    ps == vs && return ps, vs
    if ps > vs
        return max(vs, (ps ÷ vs) * vs), vs
    else
        return ps, max(ps, (vs ÷ ps) * ps)
    end
end

function _auto_amr_align_cached_intervals!(ka::KA, ps_interval, vs_interval)
    _auto_amr_should_align_intervals(ka, ps_interval, vs_interval) || return nothing
    status = ka.kinfo.status
    ps, vs = _auto_amr_multiple_pair(status.ps_interval_cached, status.vs_interval_cached)
    status.ps_interval_cached = ps
    status.vs_interval_cached = vs
    return nothing
end

@inline function _should_adapt(step::Integer, interval::Integer, converge_ratio::Integer)
    return step > interval * converge_ratio
end

function _should_adapt_pair(ka::KA, ps_interval::Integer, vs_interval::Integer,
                            converge_ratio::Integer, ps_auto, vs_auto)
    ps_due = ka.kinfo.config.solver.PS_DYNAMIC_AMR &&
             _should_adapt(ka.kinfo.status.ps_adapt_step, ps_interval, converge_ratio)
    vs_due = ka.kinfo.config.solver.VS_DYNAMIC_AMR &&
             _should_adapt(ka.kinfo.status.vs_adapt_step, vs_interval, converge_ratio)
    if _auto_amr_should_align_intervals(ka, ps_auto, vs_auto)
        # In aligned-auto mode, a short-period AMR pass waits for the matching long-period pass.
        # This preserves separate PS/VS interval estimates while ensuring recovery is paid once
        # for the joint event instead of once for each staggered event.
        joint_due = ps_due && vs_due
        return joint_due, joint_due
    end
    return ps_due, vs_due
end

@inline _partition_work_weight(::InsideSolidData) = 0.0
@inline function _partition_work_weight(ps_data::PsData)
    vs_num = ps_data.vs_data.vs_num
    ps_data.bound_enc < 0 && return 2.0 * vs_num
    if ps_data.bound_enc > 0
        solid_faces = count(face -> !isempty(face) && isa(face[1], SolidNeighbor),
                            ps_data.neighbor.data)
        return (1.0 + solid_faces) * vs_num
    end
    return float(vs_num)
end
@inline _partition_work_weight(ps_data) = 0.0

"""
$(TYPEDSIGNATURES)
Return the MPI-wide load imbalance used by the partition threshold.

The metric is `max(local_weight) / mean(local_weight) - 1`, where the local weight mirrors
[`partition_weight`](@ref): velocity-cell count for fluid cells, double weight for solid cells,
and extra weight for cut cells adjacent to immersed-boundary solid neighbors.  A value of `0.10`
therefore means the heaviest rank owns about 10% more work than the mean rank.
"""
function partition_load_imbalance(ka::KA)
    local_weight = 0.0
    for tree in ka.kdata.field.trees.data
        for ps_data in tree
            local_weight += _partition_work_weight(ps_data)
        end
    end
    total_weight = MPI.Allreduce(local_weight, MPI.SUM, MPI.COMM_WORLD)
    max_weight = MPI.Allreduce(local_weight, MPI.MAX, MPI.COMM_WORLD)
    mean_weight = total_weight / MPI.Comm_size(MPI.COMM_WORLD)
    mean_weight > 0.0 || return 0.0
    return max_weight / mean_weight - 1.0
end

function _should_partition(ka::KA, interval::Integer, converge_ratio::Integer)
    _should_adapt(ka.kinfo.status.partition_step, interval, converge_ratio) || return false
    threshold = ka.kinfo.config.solver.PARTITION_IMBALANCE_THRESHOLD
    threshold <= 0.0 && return true
    return partition_load_imbalance(ka) > threshold
end

@inline function _auto_amr_target_ds(kinfo::KInfo{DIM}) where {DIM}
    level = kinfo.config.solver.AMR_DYNAMIC_PS_MAXLEVEL
    return ntuple(d -> (kinfo.config.geometry[2d] - kinfo.config.geometry[2d - 1]) /
                      kinfo.config.trees_num[d] / 2.0^level, DIM)
end

@inline function _auto_amr_impact_floor(kinfo::KInfo)
    # Use the same local relative mass/energy/heat-flux scale as VS-AMR.  The LSR floor is
    # deliberately smaller than the legacy contribution threshold, so influential tails can enter
    # the weighted statistic without letting numerically empty far tails set the transport speed.
    return min(kinfo.config.solver.ADAPT_COEFFI_VS_LOCAL,
               kinfo.config.solver.ADAPT_COEFFI_VS_LSR_FLOOR)
end

@inline function _auto_amr_ps_sensor_threshold(kinfo::KInfo)
    return kinfo.config.solver.AUTO_AMR_PS_SENSOR_FRACTION *
           kinfo.config.solver.ADAPT_COEFFI_PS
end

@inline _auto_amr_active_ps_cell(ps_data, kinfo) = false
@inline function _auto_amr_active_ps_cell(ps_data::PsData, kinfo::KInfo)
    ps_data.bound_enc < 0 && return false
    threshold = _auto_amr_ps_sensor_threshold(kinfo)
    threshold <= 0.0 && return true
    return ps_sensor(ps_data) > threshold
end

@inline function _auto_amr_transport_rate(midpoint::AbstractVector, target_ds)
    rate = 0.0
    @inbounds for d in eachindex(target_ds)
        rate += abs(midpoint[d]) / target_ds[d]
    end
    return rate
end

@inline function _auto_amr_log_rate_bin(rate::Real, rate_max::Real)
    rate > 0.0 || return 0
    rate_max > 0.0 || return 0
    scaled = min(float(rate) / float(rate_max), 1.0)
    bin = ceil(Int, AUTO_AMR_RATE_HIST_BINS +
                    AUTO_AMR_RATE_BINS_PER_OCTAVE * log2(scaled))
    return clamp(bin, 1, AUTO_AMR_RATE_HIST_BINS)
end

@inline function _auto_amr_log_rate_upper_edge(bin::Integer, rate_max::Real)
    bin >= AUTO_AMR_RATE_HIST_BINS && return float(rate_max)
    exponent = (bin - AUTO_AMR_RATE_HIST_BINS) / AUTO_AMR_RATE_BINS_PER_OCTAVE
    return float(rate_max) * 2.0^exponent
end

@inline function _auto_amr_add_rate_sample!(hist::AbstractVector, rate::Real, weight::Real,
                                            rate_max::Real)
    weight > 0.0 || return nothing
    bin = _auto_amr_log_rate_bin(rate, rate_max)
    bin == 0 && return nothing
    @inbounds hist[bin] += weight
    return nothing
end

function _auto_amr_quantile_rate(hist::AbstractVector, rate_max::Real, quantile::Real)
    total_weight = sum(hist)
    total_weight > 0.0 || return 0.0
    target = clamp(quantile, 0.0, 1.0) * total_weight
    cumulative = 0.0
    @inbounds for bin in eachindex(hist)
        cumulative += hist[bin]
        if cumulative >= target
            return _auto_amr_log_rate_upper_edge(bin, rate_max)
        end
    end
    return float(rate_max)
end

function _auto_amr_cell_max_rate(ps_data::PsData{DIM,NDF}, target_ds,
                                 kinfo::KInfo{DIM,NDF}) where {DIM,NDF}
    vs_data = ps_data.vs_data
    U = @view(ps_data.prim[2:1 + DIM])
    impact_floor = _auto_amr_impact_floor(kinfo)
    rate = 0.0
    @inbounds for c in 1:vs_data.vs_num
        midpoint = @view(vs_data.midpoint[c, :])
        df = @view(vs_data.df[c, :])
        ratio = local_contribution_ratio(ps_data.w, U, midpoint, df, vs_data.weight[c], kinfo)
        ratio > impact_floor ||
            continue
        rate = max(rate, _auto_amr_transport_rate(midpoint, target_ds))
    end
    return rate
end

function _auto_amr_cell_rate_histogram!(hist::AbstractVector, ps_data::PsData{DIM,NDF},
                                        target_ds, rate_max::Real,
                                        kinfo::KInfo{DIM,NDF}) where {DIM,NDF}
    vs_data = ps_data.vs_data
    U = @view(ps_data.prim[2:1 + DIM])
    impact_floor = _auto_amr_impact_floor(kinfo)
    @inbounds for c in 1:vs_data.vs_num
        midpoint = @view(vs_data.midpoint[c, :])
        df = @view(vs_data.df[c, :])
        ratio = local_contribution_ratio(ps_data.w, U, midpoint, df, vs_data.weight[c], kinfo)
        ratio > impact_floor ||
            continue
        _auto_amr_add_rate_sample!(
            hist,
            _auto_amr_transport_rate(midpoint, target_ds),
            ratio,
            rate_max,
        )
    end
    return hist
end

@inline _auto_amr_ps_max_rate(ps_data, target_ds, kinfo) = 0.0
@inline function _auto_amr_ps_max_rate(ps_data::PsData, target_ds, kinfo)
    _auto_amr_active_ps_cell(ps_data, kinfo) || return 0.0
    return max(
        _auto_amr_cell_max_rate(ps_data, target_ds, kinfo),
        _auto_amr_boundary_max_rate(ps_data, target_ds, kinfo),
    )
end

@inline _auto_amr_ps_rate_histogram!(hist, ps_data, target_ds, rate_max, kinfo) = hist
@inline function _auto_amr_ps_rate_histogram!(hist::AbstractVector, ps_data::PsData, target_ds,
                                             rate_max::Real, kinfo::KInfo)
    _auto_amr_active_ps_cell(ps_data, kinfo) || return hist
    _auto_amr_cell_rate_histogram!(hist, ps_data, target_ds, rate_max, kinfo)
    _auto_amr_boundary_rate_histogram!(hist, ps_data, target_ds, rate_max, kinfo)
    return hist
end

@inline function _auto_amr_touches_domain_face(ps_data::PsData, kinfo::KInfo, faceid::Integer)
    faceid > length(kinfo.config.geometry) && return false
    dir = div(faceid - 1, 2) + 1
    side = isodd(faceid) ? -0.5 : 0.5
    face_x = ps_data.midpoint[dir] + side * ps_data.ds[dir]
    scale = max(1.0, abs(kinfo.config.geometry[faceid]), ps_data.ds[dir])
    return abs(face_x - kinfo.config.geometry[faceid]) <= 1e-10 * scale
end

@inline function _auto_amr_discrete_maxwell!(df::AbstractVector, midpoint::AbstractVector,
                                             prim::AbstractVector, kinfo::KInfo{DIM,1}) where {DIM}
    df[1] = discrete_maxwell(midpoint, prim, kinfo)
    return df
end
@inline function _auto_amr_discrete_maxwell!(df::AbstractVector, midpoint::AbstractVector,
                                             prim::AbstractVector, kinfo::KInfo{DIM,NDF}) where {DIM,NDF}
    value = discrete_maxwell(midpoint, prim, kinfo)
    @inbounds for k in 1:NDF
        df[k] = value[k]
    end
    return df
end

function _auto_amr_boundary_max_rate(ps_data::PsData{DIM,NDF}, target_ds,
                                     kinfo::KInfo{DIM,NDF}) where {DIM,NDF}
    face_midpoint = similar(ps_data.midpoint)
    bc_df = Vector{Float64}(undef, NDF)
    nsamples = length(PS_BOUNDARY_SENSOR_FRACTIONS)^(DIM - 1)
    impact_floor = _auto_amr_impact_floor(kinfo)
    vs_data = ps_data.vs_data
    rate = 0.0

    for domain in kinfo.config.domain
        isdefined(domain, :bc) || continue
        faceid = domain.id
        _auto_amr_touches_domain_face(ps_data, kinfo, faceid) || continue
        dir = div(faceid - 1, 2) + 1
        for sample in 1:nsamples
            _ps_boundary_sample_point!(face_midpoint, ps_data, kinfo, faceid, dir, sample)
            bc_prim = _ps_domain_bc_prim(domain, face_midpoint, kinfo)
            bc_prim === nothing && continue
            length(bc_prim) == length(ps_data.prim) ||
                error("Domain boundary $faceid returned a primitive vector of length $(length(bc_prim)); expected $(length(ps_data.prim)).")
            bc_w = get_conserved(bc_prim, kinfo)
            bc_U = @view(bc_prim[2:1 + DIM])
            @inbounds for c in 1:vs_data.vs_num
                midpoint = @view(vs_data.midpoint[c, :])
                df = _auto_amr_discrete_maxwell!(bc_df, midpoint, bc_prim, kinfo)
                ratio = local_contribution_ratio(bc_w, bc_U, midpoint, df, vs_data.weight[c], kinfo)
                ratio > impact_floor ||
                    continue
                rate = max(rate, _auto_amr_transport_rate(midpoint, target_ds))
            end
        end
    end
    return rate
end

function _auto_amr_boundary_rate_histogram!(hist::AbstractVector, ps_data::PsData{DIM,NDF},
                                            target_ds, rate_max::Real,
                                            kinfo::KInfo{DIM,NDF}) where {DIM,NDF}
    face_midpoint = similar(ps_data.midpoint)
    bc_df = Vector{Float64}(undef, NDF)
    nsamples = length(PS_BOUNDARY_SENSOR_FRACTIONS)^(DIM - 1)
    impact_floor = _auto_amr_impact_floor(kinfo)
    vs_data = ps_data.vs_data

    for domain in kinfo.config.domain
        isdefined(domain, :bc) || continue
        faceid = domain.id
        _auto_amr_touches_domain_face(ps_data, kinfo, faceid) || continue
        dir = div(faceid - 1, 2) + 1
        for sample in 1:nsamples
            _ps_boundary_sample_point!(face_midpoint, ps_data, kinfo, faceid, dir, sample)
            bc_prim = _ps_domain_bc_prim(domain, face_midpoint, kinfo)
            bc_prim === nothing && continue
            length(bc_prim) == length(ps_data.prim) ||
                error("Domain boundary $faceid returned a primitive vector of length $(length(bc_prim)); expected $(length(ps_data.prim)).")
            bc_w = get_conserved(bc_prim, kinfo)
            bc_U = @view(bc_prim[2:1 + DIM])
            @inbounds for c in 1:vs_data.vs_num
                midpoint = @view(vs_data.midpoint[c, :])
                df = _auto_amr_discrete_maxwell!(bc_df, midpoint, bc_prim, kinfo)
                ratio = local_contribution_ratio(bc_w, bc_U, midpoint, df, vs_data.weight[c], kinfo)
                ratio > impact_floor ||
                    continue
                # Boundary states are sampled at several points along a face.  Divide by the
                # sample count so a boundary face contributes as one representative face state,
                # not as `nsamples` independent copies of the same macroscopic influence.
                _auto_amr_add_rate_sample!(
                    hist,
                    _auto_amr_transport_rate(midpoint, target_ds),
                    ratio / nsamples,
                    rate_max,
                )
            end
        end
    end
    return hist
end

"""
$(TYPEDSIGNATURES)
Estimate the relevant physical-space transport rate carried by the current velocity grids. This
is the statistic used by automatic AMR intervals.

The estimate deliberately does **not** use the single fastest velocity cell in the whole domain.
It first restricts the physical cells to the current PS-AMR active band,
`ps_sensor > AUTO_AMR_PS_SENSOR_FRACTION * ADAPT_COEFFI_PS`; smooth cells are ignored because their
fast velocity tails do not need to schedule near-term physical refinement.  Inside those active
cells, each velocity node `ξ` contributes a rate `sum(abs(ξ[d]) / Δx_target[d])`, weighted by its
relative contribution to local mass/internal-energy/heat-flux moments.  The returned rate is the
MPI-wide weighted `AUTO_AMR_RATE_QUANTILE` of these samples, so an isolated very fast tail cell
cannot shorten the interval unless enough macroscopic contribution sits in that high-speed tail.

Domain boundary states with a supplied `bc` are sampled by evaluating their Maxwellian on adjacent
active cells' velocity grids, so inlet waves can still shorten the next interval before they enter
the interior solution.
"""
function kinetic_amr_transport_rate(ka::KA{DIM,NDF}) where {DIM,NDF}
    target_ds = _auto_amr_target_ds(ka.kinfo)
    local_max_rate = 0.0
    for tree in ka.kdata.field.trees.data
        for ps_data in tree
            local_max_rate = max(local_max_rate,
                                 _auto_amr_ps_max_rate(ps_data, target_ds, ka.kinfo))
        end
    end
    rate_max = MPI.Allreduce(local_max_rate, MPI.MAX, MPI.COMM_WORLD)
    rate_max > 0.0 || return 0.0

    local_hist = zeros(Float64, AUTO_AMR_RATE_HIST_BINS)
    for tree in ka.kdata.field.trees.data
        for ps_data in tree
            _auto_amr_ps_rate_histogram!(local_hist, ps_data, target_ds, rate_max, ka.kinfo)
        end
    end
    global_hist = MPI.Allreduce(local_hist, MPI.SUM, MPI.COMM_WORLD)
    return _auto_amr_quantile_rate(
        global_hist,
        rate_max,
        ka.kinfo.config.solver.AUTO_AMR_RATE_QUANTILE,
    )
end

function _auto_interval_from_rate(rate::Real, Δt::Real, travel_fraction::Real, max_interval::Integer)
    if !(isfinite(rate) && isfinite(Δt)) || rate <= 0.0 || Δt <= 0.0
        return Int(max_interval)
    end
    interval = floor(Int, travel_fraction / (rate * Δt))
    return clamp(interval, AUTO_AMR_MIN_INTERVAL, Int(max_interval))
end

"""
$(TYPEDSIGNATURES)
Refresh cached automatic AMR intervals in `ka.kinfo.status`.

This is a lightweight scheduling helper for [`adaptive_mesh_refinement!`](@ref). If neither
`ps_interval` nor `vs_interval` is `:auto`, it returns without doing the kinetic statistic.
Otherwise, it calls [`kinetic_amr_transport_rate`](@ref), stores the rate in
`status.amr_transport_rate`, and updates the requested cached intervals. The update is normally
performed after AMR/recovery, so the MPI-wide statistic is not part of the per-step path.

The cached intervals start from a conservative shared setup value (`5` for both physical and
velocity space; see [`Status`](@ref)) and are then estimated as
`floor(travel_fraction / (rate * Δt))`. The travel fractions and maximum clamps are solver
configuration fields (`AUTO_AMR_PS_TRAVEL_FRACTION`, `AUTO_AMR_VS_TRAVEL_FRACTION`,
`AUTO_AMR_PS_MAX_INTERVAL`, `AUTO_AMR_VS_MAX_INTERVAL`). The statistic itself is controlled by
`AUTO_AMR_PS_SENSOR_FRACTION` and `AUTO_AMR_RATE_QUANTILE`, so small calibration sweeps can tune
the accuracy/cost trade-off without editing this file.

When both PS and VS intervals are `:auto` and `AUTO_AMR_ALIGN_INTERVALS=true`, the two estimates
are adjusted to an integer-multiple pair.  The shorter estimate is kept, and the longer estimate
is rounded down to the nearest multiple of the shorter one.  This never relaxes the
transport-based interval estimate, while giving the scheduler a clean pair of cadences that can
be co-triggered without standalone [`amr_recover!`](@ref) calls.
"""
function refresh_auto_amr_intervals!(ka::KA; ps_interval = :auto, vs_interval = :auto)
    ps_auto = ps_interval === :auto
    vs_auto = vs_interval === :auto
    (ps_auto || vs_auto) || return nothing
    status = ka.kinfo.status
    solver = ka.kinfo.config.solver
    rate = kinetic_amr_transport_rate(ka)
    status.amr_transport_rate = rate
    Δt = status.Δt_ξ
    if ps_auto
        status.ps_interval_cached =
            _auto_interval_from_rate(rate, Δt, solver.AUTO_AMR_PS_TRAVEL_FRACTION,
                                     solver.AUTO_AMR_PS_MAX_INTERVAL)
    end
    if vs_auto
        status.vs_interval_cached =
            _auto_interval_from_rate(rate, Δt, solver.AUTO_AMR_VS_TRAVEL_FRACTION,
                                     solver.AUTO_AMR_VS_MAX_INTERVAL)
    end
    _auto_amr_align_cached_intervals!(ka, ps_interval, vs_interval)
    return nothing
end

"""
$(TYPEDSIGNATURES)
Run one interval-driven AMR scheduling pass.

This is the AMR stage used by [`solve!`](@ref). It may run physical-space AMR, velocity-space AMR,
partitioning, and the required recovery depending on the counters stored in `ka.kinfo.status`.
If no interval is due, the function returns without changing the mesh.

# Interval keywords

- `ps_interval=:auto`: physical-space AMR interval.
- `vs_interval=:auto`: velocity-space AMR interval.
- `partition_interval=:auto`: load-balancing interval. The default means twice the currently
  resolved `ps_interval` (for example, `ps_interval=20` gives `partition_interval=40`). A positive
  integer fixes the interval explicitly.

For `ps_interval` and `vs_interval`, accepted forms are:

- `:auto` (default): use cached automatic intervals. The cache starts short and is refreshed only
  after AMR/recovery using [`kinetic_amr_transport_rate`](@ref), avoiding an MPI-wide statistic at
  every step.
  The travel fractions and maximum clamps are configured on [`Solver`](@ref) with
  `AUTO_AMR_PS_TRAVEL_FRACTION`, `AUTO_AMR_VS_TRAVEL_FRACTION`,
  `AUTO_AMR_PS_MAX_INTERVAL`, and `AUTO_AMR_VS_MAX_INTERVAL`.
  The kinetic transport statistic is additionally gated by
  `AUTO_AMR_PS_SENSOR_FRACTION * ADAPT_COEFFI_PS` and reduced with the contribution-weighted
  `AUTO_AMR_RATE_QUANTILE`, so smooth-region high-speed tails and isolated fast velocity cells do
  not force unnecessary AMR/recovery passes.
  If both PS and VS intervals are automatic, `AUTO_AMR_ALIGN_INTERVALS=true` adjusts them to an
  integer-multiple pair and waits until both are due before running either AMR stage.  This avoids
  standalone recovery passes while preserving separate PS/VS interval estimates.  Set it to
  `false` to recover independent PS/VS automatic intervals.
- Positive integer: fixed number of steps between AMR checks, matching the legacy behavior.
- Function: custom interval policy. Accepted call signatures are `(p4est, ka, kind)`,
  `(p4est, ka)`, `(ka, kind)`, or `(ka)`, where `kind` is `:ps` or `:vs`. The function may return
  a positive integer or `:auto`.

Partitioning is checked only after a physical-space or velocity-space AMR pass has actually run.
If no AMR occurred in the current scheduler call, no load balancing is attempted even when
`partition_step` is larger than the resolved partition interval.  A due partition is also gated by
`PARTITION_IMBALANCE_THRESHOLD`: the weighted load imbalance must exceed this threshold before
`ps_partition!` is called.  Set the threshold to `0` to recover the old interval-only behavior.

The convergence safeguard used by the previous scheduler is preserved: the effective interval is
multiplied by a residual-dependent `converge_ratio` near convergence.

# Other keywords

- `ps_recursive=false`: pass recursive refinement/coarsening to physical-space AMR.
- `vs_balance=false`: balance neighboring velocity-space grids during VS-AMR.
"""
function adaptive_mesh_refinement!(p4est::P_pxest_t,ka::KA;ps_interval=:auto,vs_interval=:auto,partition_interval=:auto,ps_recursive = false, vs_balance=false)
    ka.kinfo.status.residual.redundant_step>0&&(return nothing)
    res = maximum(ka.kinfo.status.residual.residual)
    converge_ratio = res/ka.kinfo.config.solver.TOLERANCE>100 ? 1 : Int(floor(100*ka.kinfo.config.solver.TOLERANCE/res))
    _auto_amr_align_cached_intervals!(ka, ps_interval, vs_interval)
    ps_interval_value = _interval_value(ps_interval, p4est, ka, :ps)
    vs_interval_value = _interval_value(vs_interval, p4est, ka, :vs)
    partition_interval_value = _partition_interval_value(partition_interval, ps_interval_value)
    ps_changed = false
    vs_changed = false
    partition_changed = false
    ps_due, vs_due =
        _should_adapt_pair(ka, ps_interval_value, vs_interval_value, converge_ratio,
                           ps_interval, vs_interval)
    if ps_due
        ps_adaptive_mesh_refinement!(p4est,ka;recursive = ps_recursive)
        ps_changed = true;ka.kinfo.status.ps_adapt_step = 1
    end
    if vs_due
        if vs_balance && ps_changed
            update_ghost!(p4est,ka)
            update_neighbor!(p4est,ka)
        end
        vs_changed = vs_adaptive_mesh_refinement!(ka;vs_balance = vs_balance)
        ka.kinfo.status.vs_adapt_step=1
    end
    if (ps_changed || vs_changed) &&
        _should_partition(ka, partition_interval_value, converge_ratio)
        ps_partition!(p4est, ka)
        partition_changed = true;ka.kinfo.status.partition_step = 1
    end
    if ps_changed || partition_changed
        amr_recover!(p4est,ka;topology_changed = true, velocity_changed = vs_changed)
    elseif vs_changed
        amr_recover!(p4est,ka;topology_changed = false, velocity_changed = true)
    end
    # Automatic intervals are intentionally updated only after adaptation/recovery.  This keeps
    # the all-rank kinetic statistic off the per-step path: `:auto` begins with the short cached
    # interval from `Status`, then the refreshed transport estimate can safely relax it.
    (ps_changed || vs_changed || partition_changed) &&
        refresh_auto_amr_intervals!(ka; ps_interval = ps_interval, vs_interval = vs_interval)
    return nothing
end
function update_faces!(p4est::P_pxest_t, ka::KA)
    empty!(ka.kdata.field.faces)
    initialize_faces!(p4est, ka)
end

_has_immersed_boundaries(ka::KA) = !isempty(ka.kinfo.config.IB)

function _finish_immersed_boundary_recover!(ka::KA)
    if _has_immersed_boundaries(ka)
        initialize_immersed_boundaries!(ka)
    else
        empty!(ka.kdata.field.immersed_boundaries)
    end
    return nothing
end

"""
$(TYPEDSIGNATURES)
Recover the ghost layers, neighbor relations, immersed boundaries, and faces after an AMR or partition process.
"""
function amr_recover!(p4est::P_pxest_t,ka::KA;topology_changed::Bool = true, velocity_changed::Bool = true)
    if topology_changed
        update_ghost!(p4est, ka)
        update_neighbor!(p4est, ka)
    elseif velocity_changed && MPI.Comm_size(MPI.COMM_WORLD) > 1
        vs_ghost_exchange!(p4est, ka)
    end

    has_ib = _has_immersed_boundaries(ka)
    has_ib && update_solid!(ka)

    if topology_changed || has_ib || (velocity_changed && MPI.Comm_size(MPI.COMM_WORLD) > 1)
        update_faces!(p4est, ka)
    end
    _finish_immersed_boundary_recover!(ka)
    return nothing
end
