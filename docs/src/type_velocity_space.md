# Velocity space
In kinetic methods, variables are also discretized in velocity space. Information in velocity space are separately stored in [`AbstractPsData`](@ref), and are organized as following.

!!! note "Origin must be on a grid corner"
    The velocity grid is set up from the `quadrature` range and `vs_trees_num` of the
    [`Configure`](@ref). It carries an implicit requirement that the origin `v = 0` lie on a grid
    **corner** (never inside a cell); see [`check_vs_setting`](@ref) and the note on the
    [Configuration](@ref) page. The requirement is validated when the configuration is built.
```@docs
AbstractVsData
```
```@docs
VsData
```
```@docs
GhostVsData
```
```@docs
FaceVsData
```
```@docs
CuttedVelocityCells
```