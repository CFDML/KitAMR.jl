module KitAMR

using MPI
using JLD2
using Dates
using LinearAlgebra
using Parameters
using SpecialFunctions
using StaticArrays
# using PythonCall
using Suppressor
using WriteVTK
using Documenter
using Reexport
using CSV
using DataFrames

include("../lib/P4est/src/P4est.jl")
@reexport using .P4est
include("../lib/P4estTypes/src/P4estTypes.jl")
using .P4estTypes
include("../lib/KitCore/KitCore.jl")
using .KitCore
include("abstract.jl")
include("auxiliary.jl")
include("dim.jl")
include("./gas/types.jl")
include("./velocity_space/types.jl")
include("physical_space/types.jl")
include("boundary/types.jl")
include("solver/types.jl")
include("./p4est/p4est_wrap.jl")
include("model.jl")
include("math.jl")
include("IO/types.jl")
include("IO/IO.jl")
include("connectivity.jl")
include("boundary/boundary.jl")
include("adaptive.jl")
include("neighbor.jl")
include("ghost.jl")
include("velocity_space/vs_space.jl")
include("velocity_space/vs_adaptive.jl")
include("initialize.jl")
include("partition.jl")
include("slope.jl")
# include("boundary/maxwellian.jl")
include("flux/flux.jl")
include("iterate.jl")
include("finalize.jl")

# const np = Ref{Py}()
# const scipy = Ref{Py}()
# function __init__()
#     np[] = pyimport("numpy")
#     scipy[] = pyimport("scipy")
# end

end # module
