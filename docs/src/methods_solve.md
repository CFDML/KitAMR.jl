# Running a simulation

The whole lifecycle of a run is four calls:

```julia
p4est, ka = initialize(config)   # build the solver state from a `Configure`
solve!(p4est, ka)                # run the time-stepping loop to completion
save_result(p4est, ka)           # write the converged field
finalize!(p4est, ka)             # release the C-managed memory
```

[`initialize`](@ref) is documented under [Initialization](@ref), [`save_result`](@ref) under
[IO](@ref), and [`finalize!`](@ref) under [Finalization](@ref). The callbacks referenced from a
[`Configure`](@ref) (initial condition, refinement flags, velocity-space output) are collected on
the [User-defined functions](@ref) page.

`solve!` encapsulates the per-step loop and exposes all of its parameters as keyword arguments:

```@docs
solve!
```

Each stage that `solve!` runs is itself a public function, so you can write a custom loop instead
of (or around) `solve!`: [`adaptive_mesh_refinement!`](@ref), [`limit_Δt!`](@ref),
[`slope!`](@ref), [`flux!`](@ref), [`iterate!`](@ref), [`check_for_animsave!`](@ref),
[`check!`](@ref), [`check_for_convergence`](@ref), [`reached_max_time`](@ref).
