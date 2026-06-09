# Configuration

To simulate a problem using KitAMR, user is required to define the problem. Most of the user-provided information is collected in `Configure` struct. 

```@docs
Configure
```

!!! note "Velocity-space origin must be on a grid corner"
    The Cartesian velocity space carries an implicit requirement: the origin `v = 0` must lie on a
    velocity-grid **corner** (a root-cell vertex), never inside a cell. This keeps the upwinding
    split at `v·n = 0` unambiguous, and since refinement only bisects cells, 0 then stays on a
    corner at every level. Concretely, for each dimension `(0 - min)/(max - min)*vs_trees_num` must
    be an integer — easiest with a symmetric range `[-a, a]` and an even `vs_trees_num`. Satisfying
    it is left to the user (it is trivial); the setting is validated when a [`Configure`](@ref) is
    constructed, and an informative error is thrown if it is violated.

```@docs
check_vs_setting
```

```@docs
Solver
```

```@docs
UDF
```