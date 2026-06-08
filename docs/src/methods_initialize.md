# Initialization
```@docs
initialize
```
The initial condition can be re-applied to a (re)refined mesh — e.g. between the initial
pre-refinement passes — with
```@docs
reinitialize_initial_condition!
```
---
```@docs
initialize_ghost
```
---
```@docs
initialize_trees!
```
In this function, `p4est` is initialized. Mesh (in both physical and velocity space) is generated and adaptively refined. Field (in both physical and velocity space) is initialized according to initial conditions. Mapping between `Julia` data and `p4est` is established.
Currently, realizations for 2D and 3D cases are different. The method for 3D case can be more complicated but more efficient. The method for 2D will be improved to be consistent with the 3D case.

---
```@docs
initialize_ps!
```
Initialize [`PsData`](@ref) and establish the mapping between `Julia` data and `p4est`.

---
```@docs
pre_refine!
```
Mesh identification, geometric adaptive refinement, and initial partition.

---
```@docs
initialize_faces!
```
Initialize faces defined by [`AbstractFace`](@ref).

---
```@docs
AMR_ghost_new
```
Initialize the ghost layer in `p4est`.

---

```@docs
AMR_mesh_new
```
Initialize the mesh in `p4est` that is used to look up the neighboring cell in physical space.

---

```@docs
initialize_neighbor_data!
```
Initialize [`Neighbor`](@ref) in [`PsData`](@ref).

---
```@docs
initialize_solid_neighbor!
```
Initialize [`SolidNeighbor`](@ref)s corresponding to the immersed boundaries.
