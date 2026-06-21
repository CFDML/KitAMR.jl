# Automatic AMR intervals use a short cached startup value from `Status`; after an AMR pass
# this file refreshes the cache from scheduler diagnostics and may relax it.
const AUTO_AMR_MIN_INTERVAL = 1
const AUTO_AMR_MAX_INTERVAL = 50
const AUTO_AMR_PS_SENSOR_FRACTION = 0.5
const AUTO_AMR_RATE_QUANTILE = 0.8
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

@inline function _should_adapt(step::Integer, interval::Integer, converge_ratio::Integer)
    return step > interval * converge_ratio
end

function _should_adapt_pair(ka::KA, ps_interval::Integer, vs_interval::Integer,
                            converge_ratio::Integer)
    ps_due = ka.kinfo.config.solver.AMR_PS_DYNAMIC &&
             _should_adapt(ka.kinfo.status.ps_adapt_step, ps_interval, converge_ratio)
    vs_due = ka.kinfo.config.solver.AMR_VS_DYNAMIC &&
             _should_adapt(ka.kinfo.status.vs_adapt_step, vs_interval, converge_ratio)
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
    return partition_load_imbalance(ka) > PARTITION_IMBALANCE_THRESHOLD
end

@inline function _auto_amr_target_ds(kinfo::KInfo{DIM}) where {DIM}
    level = kinfo.config.solver.AMR_PS_DYNAMIC_MAXLEVEL
    return ntuple(d -> (kinfo.config.geometry[2d] - kinfo.config.geometry[2d - 1]) /
                      kinfo.config.trees_num[d] / 2.0^level, DIM)
end

@inline function _auto_amr_impact_floor(kinfo::KInfo)
    # Use the same local relative mass/energy/heat-flux floor as the LSR VS-AMR criterion.
    return kinfo.config.solver.AMR_VS_CONTRI_FLOOR
end

@inline function _auto_amr_ps_sensor_threshold(kinfo::KInfo)
    return AUTO_AMR_PS_SENSOR_FRACTION * kinfo.config.solver.AMR_PS_THRES
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
        rate = max(rate, abs(midpoint[d]) / target_ds[d])
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

@inline _auto_amr_ps_quantile_rate!(hist, ps_data, target_ds, rate_max, quantile, kinfo) = 0.0
function _auto_amr_ps_quantile_rate!(hist::AbstractVector, ps_data::PsData, target_ds,
                                     rate_max::Real, quantile::Real, kinfo::KInfo)
    _auto_amr_active_ps_cell(ps_data, kinfo) || return 0.0
    fill!(hist, 0.0)
    _auto_amr_ps_rate_histogram!(hist, ps_data, target_ds, rate_max, kinfo)
    return _auto_amr_quantile_rate(hist, rate_max, quantile)
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
is the kinetic transport statistic used by the shared automatic AMR interval.

The estimate deliberately does **not** use the single fastest velocity cell in the whole domain.
It first restricts the physical cells to the current PS-AMR active band,
`ps_sensor > 0.5 * AMR_PS_THRES`; smooth cells are ignored because their
fast velocity tails do not need to schedule near-term physical refinement.  Inside those active
cells, each velocity node `ξ` contributes a rate `max(abs(ξ[d]) / Δx_target[d])`, weighted by its
relative contribution to local mass/internal-energy/heat-flux moments.  The statistic first takes
the contribution-weighted `0.8` quantile inside each active physical cell, then returns the
MPI-wide maximum of those per-cell quantile rates.  This keeps the interval controlled by the
fastest locally relevant structure without letting the number of samples in other cells dilute it.

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

    cell_hist = zeros(Float64, AUTO_AMR_RATE_HIST_BINS)
    local_quantile_max = 0.0
    for tree in ka.kdata.field.trees.data
        for ps_data in tree
            local_quantile_max = max(
                local_quantile_max,
                _auto_amr_ps_quantile_rate!(cell_hist, ps_data, target_ds, rate_max,
                                            AUTO_AMR_RATE_QUANTILE, ka.kinfo),
            )
        end
    end
    return MPI.Allreduce(local_quantile_max, MPI.MAX, MPI.COMM_WORLD)
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
`ps_interval` nor `vs_interval` is `:auto`, it returns immediately.  Otherwise one shared cached
interval is estimated from [`kinetic_amr_transport_rate`](@ref) and copied to both PS and VS auto
caches.  This keeps automatic PS/VS scheduling on the same physical-propagation cadence.
"""
function refresh_auto_amr_intervals!(ka::KA; ps_interval = :auto, vs_interval = :auto)
    ps_auto = ps_interval === :auto
    vs_auto = vs_interval === :auto
    (ps_auto || vs_auto) || return nothing
    status = ka.kinfo.status
    solver = ka.kinfo.config.solver
    rate = kinetic_amr_transport_rate(ka)
    interval = _auto_interval_from_rate(rate, status.Δt_ξ, solver.AUTO_AMR_PS_TRAVEL_FRACTION,
                                       AUTO_AMR_MAX_INTERVAL)
    status.amr_transport_rate = rate
    status.ps_interval_cached = interval
    status.vs_interval_cached = interval
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
  after AMR/partition events or VS-AMR checks, avoiding an MPI-wide statistic at every step.
  Automatic PS-AMR and VS-AMR use the same kinetic physical-propagation estimate controlled by
  `AUTO_AMR_PS_TRAVEL_FRACTION`, with fixed internal defaults for the sensor gate (`0.5`) and
  contribution quantile (`0.8`). Automatic intervals are capped at 50 steps.
- Positive integer: fixed number of steps between AMR checks, matching the legacy behavior.
- Function: custom interval policy. Accepted call signatures are `(p4est, ka, kind)`,
  `(p4est, ka)`, `(ka, kind)`, or `(ka)`, where `kind` is `:ps` or `:vs`. The function may return
  a positive integer or `:auto`.

Partitioning is checked only after a physical-space or velocity-space AMR pass has actually run.
If no AMR occurred in the current scheduler call, no load balancing is attempted even when
`partition_step` is larger than the resolved partition interval.  A due partition is also gated by
the default weighted load-imbalance threshold before `ps_partition!` is called.

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
    ps_interval_value = _interval_value(ps_interval, p4est, ka, :ps)
    vs_interval_value = _interval_value(vs_interval, p4est, ka, :vs)
    partition_interval_value = _partition_interval_value(partition_interval, ps_interval_value)
    ps_changed = false
    vs_changed = false
    vs_checked = false
    partition_changed = false
    ps_due, vs_due =
        _should_adapt_pair(ka, ps_interval_value, vs_interval_value, converge_ratio)
    if ps_due
        ps_adaptive_mesh_refinement!(p4est,ka;recursive = ps_recursive)
        ps_changed = true;ka.kinfo.status.ps_adapt_step = 1
    end
    if vs_due
        vs_checked = true
        if vs_balance && ps_changed
            update_ghost!(p4est,ka)
            update_neighbor!(p4est,ka)
        end
        vs_changed = _vs_adaptive_mesh_refinement_result!(ka;vs_balance = vs_balance)
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
    # Automatic intervals are intentionally updated only after AMR/partition events or an actual
    # VS-AMR check.  This keeps the all-rank statistics off the per-step path.
    (ps_changed || vs_checked || partition_changed) &&
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

    if topology_changed || has_ib
        update_faces!(p4est, ka)
    end
    _finish_immersed_boundary_recover!(ka)
    return nothing
end
