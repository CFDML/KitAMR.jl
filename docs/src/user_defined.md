# User-defined functions

A simulation is configured almost entirely with data, but a few behaviours are supplied as
**callbacks** — small functions you write and hand to KitAMR through the configuration objects.
They fall into three groups:

| Callback | Configured via | Purpose |
|---|---|---|
| Initial condition | [`PCoordFn`](@ref) (in `Configure(; IC = …)`) | the primitive state at each cell centre |
| Refinement flags | [`UDF`](@ref) (in `Configure(; user_defined = …)`) | steer adaptive mesh refinement |
| Velocity-space output | [`Output`](@ref) (in `Configure(; output = …)`) | choose which cells dump their velocity space |

All callbacks are **optional**: every field defaults to a no-op, so a uniform run needs none of
them. They are collected here because they share conventions (the primitive-variable ordering,
the cell objects passed in) and because getting their signatures right is the only non-obvious
part of writing a case.

## The primitive-variable vector

Initial and boundary states are given as a **primitive vector**

```
[ρ, U₁, …, U_DIM, λ]          # length DIM + 2
```

where `ρ` is density, `Uᵢ` the bulk-velocity components, and `λ = ρ/(2p)` the inverse
temperature (`p` is the pressure). In 2D this is `[ρ, u, v, λ]` (length 4); in 3D
`[ρ, u, v, w, λ]` (length 5). The same ordering is used by [`Uniform`](@ref), [`PCoordFn`](@ref)
and the [`Domain`](@ref) boundary states.

## Initial condition

Use [`Uniform`](@ref) for a constant state, or [`PCoordFn`](@ref) for a coordinate-dependent one.
`PCoordFn` wraps a function

```julia
PCIC_fn(midpoint::Vector{Float64}, kinfo::KInfo) -> Vector{Float64}   # the primitive vector
```

- `midpoint` — cell-centre coordinate (length `DIM`).
- `kinfo` — the [`KInfo`](@ref) object (gives access to gas properties, geometry, …).
- **returns** the primitive vector `[ρ, U₁, …, U_DIM, λ]` for that cell.

It is evaluated once per cell during [`initialize`](@ref); when newly refined cells need the
exact initial state again, [`reinitialize_initial_condition!`](@ref) re-evaluates it (this is what
`solve!`'s `prerefine_reinit_ic` does).

```julia
# 2D buffer initial condition around a cylinder (from example/cylinder)
function cylinder_buffer_IC(midpoint::Vector{Float64}, ::KInfo)
    r = norm(midpoint); Ma = 5.0; Tw = 1.0; R = 1.0; l = 1.0   # buffer length
    if r > R + l
        return [1.0, Ma*√(5/6), 0., 1.0]
    else
        return [1.0, (r-R)*Ma*√(5/6), 0., Tw-(r-R)/l*(Tw-1.0)]
    end
end

config = Configure(solver; IC = PCoordFn(cylinder_buffer_IC), ...)
```

## Refinement flags ([`UDF`](@ref))

Two flags steer physical-space adaptive mesh refinement. Supply either or both via
`UDF(; static_ps_refine_flag = …, dynamic_ps_refine_flag = …)` and pass it as
`Configure(; user_defined = udf)`.

### `static_ps_refine_flag` — geometry-driven, applied once at setup

```julia
static_ps_refine_flag(midpoint::Vector{Float64}, ds::Vector{Float64},
                      kinfo::KInfo, level::Int) -> Bool
```

Evaluated during the initial refinement ([`pre_refine!`](@ref), inside [`initialize`](@ref)).
`midpoint`/`ds` are the candidate cell's centre/size and `level` its current level. Return `true`
to **force-refine** the cell; force-refined cells are also protected from later coarsening. Use it
to guarantee resolution in a fixed region (e.g. a shock layer around a body) regardless of the
solution. Default: never force.

### `dynamic_ps_refine_flag` — solution-driven, applied every step

```julia
dynamic_ps_refine_flag(ps_data::AbstractPsData, level::Int, ka::KA) -> Bool
```

Evaluated every refinement step ([`adaptive_mesh_refinement!`](@ref)) to **gate** the built-in
Löhner sensor. Return `false` to forbid dynamic refinement of `ps_data` (e.g. to freeze the mesh
outside a region of interest); `true` lets the sensor decide. The cell object exposes
`ps_data.midpoint`, `ps_data.ds`, `ps_data.prim`, … Default: always allow.

```julia
# from example/cylinder
shock_wave_region(midpoint, ds, kinfo, level) =
    (-5 < midpoint[1] < 5 && -5 < midpoint[2] < 5 && √sum(midpoint.^2) > 1.0 && level < 4)

amr_region(ps_data, level, ka) =
    (m = ps_data.midpoint; -5 < m[1] < 5 && -5 < m[2] < 5 && √sum(m.^2) > 1.0)

udf = UDF(; static_ps_refine_flag = shock_wave_region,
            dynamic_ps_refine_flag = amr_region)
config = Configure(solver; user_defined = udf, ...)
```

!!! note
    The `static_vs_refine_flag` field of [`UDF`](@ref) is reserved for a future
    velocity-space refinement hook and is **currently unused**.

## Velocity-space output ([`Output`](@ref))

By default KitAMR writes the macroscopic flow field plus, on request, the full per-cell velocity
distribution for selected cells. `output.vs_output_criterion` chooses those cells. It has **two
call conventions**, one per output path; write whichever matches the output you use:

- Final result ([`save_result`](@ref)):

  ```julia
  vs_output_criterion(ps_data, ka) -> Bool          # true ⇒ include this cell
  ```
  Default (unset): include every cell.

- Animation frames ([`check_for_animsave!`](@ref), enabled by `Output(; anim_dt = …)`):

  ```julia
  vs_output_criterion(; ps_data, ka) -> (id::Int, flag::Bool)   # flag selects, id names the file
  ```
  Default (unset): write no per-cell velocity space.

```julia
# Animation convention: track 5 probe points (from example/cylinder)
function vs_output_flag(; ps_data, kwargs...)
    midpoint = ps_data.midpoint; ds = ps_data.ds
    tps = [[-1.459, 0.057], [-1.379, 0.057], [1.731, 3.347],
           [-0.904, 0.454], [0.936, 0.979]]
    for i in eachindex(tps)
        if abs(midpoint[1]-tps[i][1]) < 0.5ds[1] && abs(midpoint[2]-tps[i][2]) < 0.5ds[2]
            return i, true
        end
    end
    return 0, false
end

output = Output(solver; vs_output_criterion = vs_output_flag)
config = Configure(solver; output = output, ...)
```

## Putting it together

```julia
udf    = UDF(; static_ps_refine_flag = shock_wave_region,
               dynamic_ps_refine_flag = amr_region)
output = Output(solver; vs_output_criterion = vs_output_flag)
config = Configure(solver;
    IC           = PCoordFn(cylinder_buffer_IC),   # initial condition callback
    user_defined = udf,                            # refinement-flag callbacks
    output       = output,                         # velocity-space output callback
    ...)
p4est, ka = initialize(config)
solve!(p4est, ka)
```
