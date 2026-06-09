# Initial and boundary conditions

To launch a simulation, initial and boundary conditions are necessary.

## Initial conditions
The types of initial conditions are abstracted as 
```@docs
AbstractInitCond
```
Currently, available options are
```@docs
Uniform
```

```@docs
PCoordFn
```

## Boundary conditions
The type of the boundary conditions is determined by 
```@docs
AbstractBoundCond
```
Currently, available options are
```@docs
Maxwellian
```

```@docs
SuperSonicInflow
```

```@docs
SuperSonicOutflow
```

```@docs
UniformOutflow
```

```@docs
InterpolatedOutflow
```

```@docs
Period
```

With these types of boundary conditions, 
```@docs
AbstractBoundary
```
can be constructed. There are two types of the boundary.
### Domain boundary
KitAMR.jl always adopts square simulation domain. At the edges of the domain, proper conditions are required to maintain the well-posedness.

```@docs
Domain
```


### Immersed boundary
KitAMR.jl adopts the immersed boundary method (IBM) to resolve boundaries with complex geometry.
To define the geometry, KitAMR.jl provides following interface.

As the special case, the cirlce in 2D and sphere in 3D are defined separately, and are abstracted as `AbstractCircle=Union{Circle,Sphere}`.

```@docs
Circle
```
```@docs
Sphere
```

For other cases, KitAMR.jl read `.csv` file for 2D, which defines the vertices coordinates of the closed boundary curve; and `.stl` file for 3D, which provides the triangular discretization of the boundary surface. The example can be find in `/example/airfoil` and `/example/X38` respectively.
```@docs
Vertices
```
The most common constructor is
```@docs
Vertices(::Type{T},file::String,solid,refine_coeffi,bc) where{T<:KitAMR.AbstractBoundCond}
```
```@docs
Triangles
```
The most common constructor is
```@docs
Triangles(::Type{T},file::String,solid,search_radius,bc) where{T<:KitAMR.AbstractBoundCond}
```
`TriangleKDT` is a struct containing information related to K-D tree for efficient mesh generation:
```@docs
TriangleKDT
```