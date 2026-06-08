# Iteration

In kinetic methods, time marching is commonly coupled with collision process. The time marching scheme as well as the realization of collision process is abstracted as
```@docs
AbstractTimeMarchingType
```

The best supported time marching type is
```@docs
CAIDVM_Marching
```
In this type, time derivative is dicretized with backward Euler method, exhibiting first order accuracy. Collision process is carried out with conserved correction, which is crucial for AMR in velocity space.

---
Here are the related functions. The per-step driver advances the solution one time step
(collision + time marching), reduces the residual across ranks, exchanges ghost data and
increments the step counters:

```@docs
iterate!(::KitAMR.P_pxest_t,::KA)
```

The time step itself is set from the CFL/grid step and, when needed, shrunk to land exactly on
the next output/termination time by
```@docs
limit_Δt!
```
and the run terminates when
```@docs
reached_max_time
```
---

Developer may hope to define their own collision process.
As an example, the method for [`CAIDVM_Marching`](@ref) is defined as
```@docs
iterate!(::Type{CAIDVM_Marching},::KA)
```