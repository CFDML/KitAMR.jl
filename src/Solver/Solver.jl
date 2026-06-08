include("AMR.jl")
include("Auxiliary.jl")
include("Finalize.jl")
include("Initialize.jl")
export Configure, Uniform, PCoordFn, Solver, Output, UDF
export KA, KInfo, KData, Forest, Status, Residual, Ghost, GhostBuffer, GhostInfo, PsTrees, Field
export residual_check!, finalize!, check_for_convergence
export initialize, initialize_ps!, initialize_ghost, initialize_trees!, pre_refine!, initialize_faces!, reinitialize_initial_condition!
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

Initial mesh adaptation (before the loop):
- `prerefine_steps::Integer=solver.AMR_DYNAMIC_PS_MAXLEVEL` — number of
  [`ps_adaptive_mesh_refinement!`](@ref) passes; the default builds the mesh up to the
  dynamic physical-space max level. Pass `0` to skip.
- `prerefine_recursive::Bool=false` — `recursive` flag for those passes.
- `prerefine_reinit_ic::Bool=true` — re-apply the initial condition
  ([`reinitialize_initial_condition!`](@ref)) after each pass, so newly refined cells get
  the exact IC (recommended for sharp initial conditions, e.g. Riemann / blast waves).

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
"""
function solve!(
        p4est::P_pxest_t, ka::KA;
        prerefine_steps::Integer = ka.kinfo.config.solver.AMR_DYNAMIC_PS_MAXLEVEL,
        prerefine_recursive::Bool = false, prerefine_reinit_ic::Bool = true,
        ps_interval = 40, vs_interval = 80, partition_interval = 40,
        ps_recursive::Bool = false, vs_balance::Bool = false,
        max_steps = typemax(Int), break_on_convergence::Bool = true,
        listen_for_save::Bool = true,
        animation::Bool = true, anim_path::String = "./animation",
        status_check::Bool = true,
    )
    for _ in 1:prerefine_steps
        ps_adaptive_mesh_refinement!(p4est, ka; recursive = prerefine_recursive)
        prerefine_reinit_ic && reinitialize_initial_condition!(ka)
        amr_recover!(p4est, ka)   # rebuild ghost/neighbor/faces after the mesh changed
    end
    listen_for_save && listen_for_save!()
    while !reached_max_time(ka)
        adaptive_mesh_refinement!(p4est, ka;
            ps_interval = ps_interval, vs_interval = vs_interval,
            partition_interval = partition_interval,
            ps_recursive = ps_recursive, vs_balance = vs_balance)
        limit_Δt!(ka)
        slope!(p4est, ka)
        flux!(p4est, ka)
        iterate!(p4est, ka)
        animation && check_for_animsave!(p4est, ka; path = anim_path)
        status_check && check!(p4est, ka)
        break_on_convergence && check_for_convergence(ka) && break
        ka.kinfo.status.step ≥ max_steps && break
    end
    return ka
end
