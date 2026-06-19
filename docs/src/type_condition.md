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
Composite
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

### Composite domain boundary

Use [`Composite`](@ref) when different parts of one Cartesian domain edge should behave like
different domain boundaries. The payload is [`CompositeBC`](@ref): it stores a weight function and
the component [`Domain`](@ref) objects to blend.

```@docs
CompositeBC
```

The weights may be constant or coordinate-dependent. They must be non-negative and are normalized
before use. During flux evaluation, each component boundary flux is computed with its usual
implementation and the results are blended. During boundary-state sampling for AMR and initial
velocity-space refinement, components without an explicit primitive boundary state (for example
[`UniformOutflow`](@ref)) are skipped.

All component domains must use the same boundary id as the outer `Domain(Composite, id, ...)`, and
each component must be supported by the selected flux scheme.

For example, an inlet aperture on the left boundary surrounded by a zero-gradient open boundary
can be written as

```julia
function aperture_weight(midpoint)
    a = clamp(0.5 * (1 - tanh((abs(midpoint[2]) - jet_radius) / edge_width)), 0.0, 1.0)
    return (a, 1.0 - a)
end

left_boundary = Domain(Composite, 1,
    CompositeBC(aperture_weight,
        Domain(SuperSonicInflow, 1, [rho_jet, u_jet, v_jet, lambda_jet]),
        Domain(UniformOutflow, 1)))
```

Here the aperture core behaves as `SuperSonicInflow`, the exterior behaves as `UniformOutflow`,
and the aperture edge is smoothly blended by `aperture_weight`.


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
