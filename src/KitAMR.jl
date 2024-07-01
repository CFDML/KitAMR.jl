module KitAMR

using MPI
using P4est
using StaticArrays
using Parameters
using SpecialFunctions
using CairoMakie
using LinearAlgebra
using HDF5
using JLD2

include("../lib/P4estTypes/src/P4estTypes.jl")
using .P4estTypes

using PythonCall
scipy = pyimport("scipy")
np = pyimport("numpy")

include("abstract.jl")
include("types.jl")
include("model.jl")
include("p4est_wrap.jl")
include("math.jl")
include("dim.jl")
include("connectivity.jl")
include("vs_space.jl")
include("adaptive.jl")
include("vs_adaptive.jl")
include("initialize.jl")
include("ghost.jl")
include("neighbor.jl")
include("reconstruct.jl")
include("flux_interp_t.jl")
include("iterate.jl")
include("partition.jl")
include("postprocess.jl")
include("finalize.jl")

end # module
