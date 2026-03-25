# Reconstruction

## Slope
In current version, the distribution function is reconstructed as a linear function. Besides the function value at center of the cell, the slope is also required. The slope is obtained by finite difference. To eliminate the oscillation near discontinuities, the vanLeer limiter is adopted. 

```@docs
vanleer
```
---

The finite difference is not that simple considering cells in templates may exhibit different refinement level, in both physical and velocity space. Here we list related methods. For more detailed explanation, one may refer to our future paper.

```@docs
diff_vs!
```
---

```@docs
update_slope!(::KA{DIM,NDF}) where{DIM,NDF}
```
Outer function to update `sdf` in [`VsData`](@ref) and `sw` in [`PsData`](@ref).

---

Here are methods corresponding to different cases in physical space. The first two argument types refer to the number of the neighboring cells in negative and positive direction along the `dir`th dimension. `Val{i}` refers to there are `i` neighboring cells across the face. `-1` represents a single cell with lower refinement level. Because the mesh is balanced, only following cases are required to consider.

```@docs
update_slope!(::Val{0},::Val{1},::PsData,::KInfo,::AbstractVector,::AbstractVector,::Integer)
update_slope!(::Val{1},::Val{0},::PsData,::KInfo,::AbstractVector,::AbstractVector,::Integer)
update_slope!(::Val{0},::Val{-1},::PsData,::KInfo,::AbstractVector,::AbstractVector,::Integer)
update_slope!(::Val{-1},::Val{0},::PsData,::KInfo,::AbstractVector,::AbstractVector,::Integer)
update_slope!(::Val{0},::KitAMR.NeighborNum,::PsData,::KInfo,::AbstractVector,::AbstractVector,::Integer)
update_slope!(::KitAMR.NeighborNum,::Val{0},::PsData,::KInfo,::AbstractVector,::AbstractVector,::Integer)
```
Update both micro and macro slopes in cells neighboring to domain edges. 

---
```@docs
update_slope!(::Val{1},::Val{1},::PsData,::KInfo,::AbstractVector,::AbstractVector,::Integer)
update_slope!(::Val{1},::Val{-1},::PsData,::KInfo,::AbstractVector,::AbstractVector,::Integer)
update_slope!(::Val{-1},::Val{1},::PsData,::KInfo,::AbstractVector,::AbstractVector,::Integer)
update_slope!(::Val{1},::KitAMR.NeighborNum,::PsData,::KInfo,::AbstractVector,::AbstractVector,::Integer)
update_slope!(::KitAMR.NeighborNum,::Val{1},::PsData,::KInfo,::AbstractVector,::AbstractVector,::Integer)
update_slope!(::Val{-1},::KitAMR.NeighborNum,::PsData,::KInfo,::AbstractVector,::AbstractVector,::Integer)
update_slope!(::KitAMR.NeighborNum,::Val{-1},::PsData,::KInfo,::AbstractVector,::AbstractVector,::Integer)
update_slope!(::KitAMR.NeighborNum,::KitAMR.NeighborNum,::PsData,::KInfo,::AbstractVector,::AbstractVector,::Integer)
update_slope!(::Val{-1},::Val{-1},::PsData,::KInfo,::AbstractVector,::AbstractVector,::Integer)
```
Update both micro and macro in internal cells.

---
```@docs
update_slope_inner_vs!
update_slope_bound_vs!
```
Update micro slope `sdf`.

---
```@docs
update_slope_inner_ps!
update_slope_bound_ps!
```
Update macro slope `sw`. Currently, `sw` is only used as an indicator for AMR in physical space. Hence, the calculation is carried out without a limiter. On the contrary, the slope with larger absolute value in two sides slopes is adopted for more sufficient refinement.

---

To preserve the positivity of the distribution function, a positivity-preserving reconstruction strategy is adopted. Detailed information is in [Flux](@ref).

## Immersed boundary method

Boundaries with complex geometry are resolved with immersed boundary method. More specifically, a sharp interface method with ghost cell is adopted. Detailed information can be found in our paper https://doi.org/10.48550/arXiv.2512.20252. KitAMR.jl provided interfaces for the reconstruction of the ghost cells
```@docs
update_solid_cell!
```
```@docs
update_solid_neighbor!
```