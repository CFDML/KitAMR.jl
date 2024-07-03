module KitAMR

using CairoMakie
using HDF5
using MPI
using JLD2
using LinearAlgebra
using Parameters
using SpecialFunctions
using StaticArrays
using PythonCall

using Reexport
@reexport using P4est

include("../lib/P4estTypes/src/P4estTypes.jl")
using .P4estTypes

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

const np = Ref{Py}()
const scipy = Ref{Py}()
function __init__()
    np[] = pyimport("numpy")
    scipy[] = pyimport("scipy")
end

end # module
