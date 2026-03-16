# Tutorial

## Installation

In your personal workstation, the installation of KitAMR.jl is quite simple. Just run
```julia
julia> ]
(//your_julia_version//) pkg> activate .
(//your_path//) pkg> add KitAMR
```
For cluster environment, preliminary procedure is required to deploy a local `p4est`, which enable KitAMR.jl to use the optimized `MPI` backend on cluster.

## Utilization

To launch a simulation, user is required to execute a `Julia` script from the command line in an `MPI` manner. Examples can be found in `/example`. Here we show one of them and provide detailed explanation.

### X-38 spacecraft
In KitAMR.jl, user-provided configuration can be input in two approaches. 

The first one is to define a `.txt` file, that is `/example/X38/configure_X38.txt`. The file will be read by [`read_config`](@ref) and returned as a dictionary. In the file, names in front of `=` are required to be consistent with the fieldnames of [`Configure`](@ref) type and the types in its field (including [`Gas`](@ref), [`Solver`](@ref), [`Output`](@ref), and [`UDF`](@ref)). The meaning of these fields can be referred to the definition of these types in `Type` pages.

The other one is to construct the [`Configure`](@ref) directly in `.jl` script.

The scripts in these two approaches are respectively given by `/example/X38/X38_text.jl` and `/example/X38/X38_julia.jl`.

The initial conditions function in [`PCoordFn`](@ref) and other user-defined functions in [`UDF`](@ref) are collected in `/example/X38/X38_udf.jl`.

The geometry of the immersed boundary is defined in `/example/X38/X38_normalized.stl`.

For other information, one may refer to the comments in `/example/X38/X38_julia.jl`.