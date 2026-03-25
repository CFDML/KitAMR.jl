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
Here are the related functions:

```@docs
iterate!(::KA)
```
---

Developer may hope to define their own collision process.
As an example, the method for [`CAIDVM_Marching`](@ref) is defined as
```@docs
iterate!(::Type{CAIDVM_Marching},::KA)
```