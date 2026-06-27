# Adaptive mesh refinement 

AMR is the core functionality of KitAMR.jl. An integrated interface is provided
```@docs
adaptive_mesh_refinement!
```

## Physical space

In physical space, a forest of quadtrees (for 2D) or octrees (for 3D) is constructed and managed with [p4est](https://github.com/cburstedde/p4est).

```@docs
ps_adaptive_mesh_refinement!
```
---
```@docs
ps_refine!
```
Refine cells satisfying criteria in physical space.

---
```@docs
ps_coarsen!
```
Coarsen cells satisfying criteria in physical space.

---
```@docs
ps_balance!
```
Balance mesh in physical space after refinement and coarsening. It ensures that the difference between the refinement levels of neighboring cells is not larger than 1.

---
In an AMR process, it is crucial to construct variables in newly generated cell (a larger one or a smaller one). Related methods are
```@docs
ps_replace!(::Val{1},out_quad,in_quads,which_tree,ka::KA{DIM}) where{DIM}
ps_replace!(::KitAMR.ChildNum,out_quads,in_quads,which_tree,ka::KA{DIM,NDF}) where{DIM,NDF}
```
---

AMR criteria in physical space are defined by
```@docs
ps_refine_flag(::PsData{DIM},level::Int8, ::KA{DIM}) where{DIM}
```
```@docs
ps_coarsen_flag(::Vector{PsData},::Vector{Int},::KA{DIM,NDF}) where{DIM,NDF}
```

By default, KitAMR.jl decides according to the physical-space Löhner sensor computed from
macroscopic variables. The dynamic criterion can be replaced with
`UDF(; dynamic_ps_adapt_criterion = ...)`; when the hook is omitted, the default Löhner criterion
is used. To compute older gradient diagnostics, a globally maximum gradient is obtained by

```@docs
update_gradmax!
```

## Velocity space

In velocity space, considering the enormous number of the phase grids, constructing a forest in velocity space for every physical cell is impossible. Hence, the AMR in velocity space is more elaborately designed.

Related functions are
```@docs
vs_adaptive_mesh_refinement!
```
```@docs
vs_refine!
```
```@docs
vs_coarsen!
```
The replacement functions are more particular and are planted inside the adaptation functions.

Compared with the AMR in physical space, the conservation is not automatically maintained in velocity space. Hence, a correction is required
```@docs
vs_conserved_correction!
```

The refinement and coarsening in velocity space use the local linear least-squares residual
indicator on the distribution function, gated by the relative local mass, energy, and heat-flux
contribution of the velocity cell.

---

## Recovery

After an AMR process, the mesh has been changed a lot. The ghost layers,  neighbor relations, immersed boundaries, and faces are required to be reconstructed.

```@docs
amr_recover!
```
Currently, the recovery function is yet to be optimized for higher efficiency. However, the overhead is negligible.
