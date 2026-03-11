# Initialization
```@docs
initialize_KitAMR
```

```@docs
initialize_ghost
```

```@docs
initialize_field!
```
In this function, `p4est` is initialized. Mesh (in both physical and velocity space) is generated and adaptively refined. Field (in both physical and velocity space) is initialized according to initial conditions. Mapping between `Julia` data and `p4est` is established.
Currently, realizations for 2D and 3D cases are different. The method for 3D case can be more complicated but more efficient. The method for 2D will be improved to be consistent with the 3D case.


Initialize [`PS_Data`](@ref) and establish the mapping between `Julia` data and `p4est`.
```@docs
initialize_ps!
```

Mesh identification, geometric adaptive refinement, and initial partition.
```@docs
pre_refine!
```

Initialize faces defined by [`AbstractFace`](@ref).
```@docs
initialize_faces!
```

Initialize the ghost layer in `p4est`.
```@docs
AMR_ghost_new
```

Initialize the mesh in `p4est` that is used to look up the neighboring cell in physical space.
```@docs
AMR_mesh_new
```

Initialize [`Neighbor`](@ref) in [`PS_Data`](@ref).
```@docs
initialize_neighbor_data!
```

Initialize [`SolidNeighbor`](@ref)s corresponding to the immersed boundaries.
```@docs
initialize_solid_neighbor!
```
