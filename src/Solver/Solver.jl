include("AMR.jl")
include("Auxiliary.jl")
include("Finalize.jl")
include("Initialize.jl")
export Configure, Uniform, PCoordFn, Solver, Output, UDF, check_vs_setting
export KA,
    KInfo, KData, Forest, Status, Residual, Ghost, GhostBuffer, GhostInfo, PsTrees, Field
export residual_check!, finalize!, check_for_convergence
export initialize,
    initialize_ps!,
    initialize_ghost,
    initialize_trees!,
    pre_refine!,
    initialize_faces!,
    reinitialize_initial_condition!
export adaptive_mesh_refinement!, amr_recover!
export solve!

"""
$(TYPEDSIGNATURES)
Run the kinetic solver's time-stepping loop to completion, so a simulation is just
`p4est, ka = initialize(config)` followed by `solve!(p4est, ka)`. Encapsulates the
canonical per-step sequence (replacing the hand-written loop in the examples):

    adaptive_mesh_refinement! → limit_Δt! → slope! → flux! → iterate!
        → check_for_animsave! → check! → (convergence break)

Returns `ka`.

# Keyword arguments

(Initial mesh pre-refinement is now done by [`initialize`](@ref) via its `prerefine_*` keyword
arguments, not here.)

AMR / load balancing (forwarded to [`adaptive_mesh_refinement!`](@ref) each step):

  - `ps_interval=40`, `vs_interval=80`, `partition_interval=40`
  - `ps_recursive::Bool=false`, `vs_balance::Bool=false`

Loop control:

  - `max_steps=typemax(Int)`       — hard cap on the number of steps.
  - `break_on_convergence::Bool=true` — stop when [`check_for_convergence`](@ref) holds.

Lifecycle / IO (toggle `false` to manage these yourself):

  - `listen_for_save::Bool=true`   — [`listen_for_save!`](@ref) before the loop.
  - `animation::Bool=true`, `anim_path::String="./animation"` — per-step [`check_for_animsave!`](@ref).
  - `status_check::Bool=true`      — per-step [`check!`](@ref) (periodic status print + save hook).
  - `progress::Bool=true`          — show a `ProgressMeter` bar (rank 0 only) whose length is one
    `ST_CHECK_INTERVAL` window, i.e. how far the current step is from the next [`check!`](@ref); it
    completes and restarts at each check.
"""
function solve!(
    p4est::P_pxest_t,
    ka::KA;
    ps_interval = 40,
    vs_interval = 80,
    partition_interval = 40,
    ps_recursive::Bool = false,
    vs_balance::Bool = false,
    max_steps = typemax(Int),
    break_on_convergence::Bool = true,
    listen_for_save::Bool = true,
    animation::Bool = true,
    anim_path::String = "./animation",
    status_check::Bool = true,
    progress::Bool = true,
)
    listen_for_save && listen_for_save!()
    # Progress bar (rank 0 only): one bar spans a single `ST_CHECK_INTERVAL` window, i.e. it shows
    # how far the current step is from the next `check!`; it completes and restarts at each check.
    interval = ka.kinfo.config.solver.ST_CHECK_INTERVAL
    bar =
        (progress && MPI.Comm_rank(MPI.COMM_WORLD) == 0) ?
        Progress(interval; desc = "Solving → next check!: ", color = :cyan) : nothing
    while !reached_max_time(ka)
        adaptive_mesh_refinement!(
            p4est,
            ka;
            ps_interval = ps_interval,
            vs_interval = vs_interval,
            partition_interval = partition_interval,
            ps_recursive = ps_recursive,
            vs_balance = vs_balance,
        )
        limit_Δt!(ka)
        slope!(p4est, ka)
        flux!(p4est, ka)
        iterate!(p4est, ka)
        animation && check_for_animsave!(p4est, ka; path = anim_path)
        atcheck = mod(ka.kinfo.status.step, interval) == 0
        # finish the bar at a window boundary so check!'s multi-line status prints cleanly below it
        if bar !== nothing && atcheck
            ProgressMeter.finish!(bar)
        end
        status_check && check!(p4est, ka)
        if bar !== nothing
            if atcheck
                bar = Progress(interval; desc = "Solving → next check!: ", color = :cyan)
            else
                ProgressMeter.update!(bar, mod(ka.kinfo.status.step, interval))
            end
        end
        break_on_convergence && check_for_convergence(ka) && break
        ka.kinfo.status.step ≥ max_steps && break
    end
    bar === nothing || ProgressMeter.finish!(bar)
    return ka
end
