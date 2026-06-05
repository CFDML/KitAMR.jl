module KitAMR

using MPI
using JLD2
using Dates
using LinearAlgebra
using Parameters
using SpecialFunctions
using StaticArrays
using Suppressor
using WriteVTK
using DocStringExtensions
using Reexport
using CSV
using DataFrames
using Statistics
using FileIO,NearestNeighbors,AbstractTrees
using GeometryBasics:Mesh
using StructArrays
using P4est

# include("../lib/P4est/src/P4est.jl")
# using .P4est
include("../lib/KitCore/KitCore.jl")
using .KitCore
include("Abstract/Types.jl")
include("Boundary/Boundary.jl")
include("Flux/Flux.jl")
include("Gas/Gas.jl")
include("IO/IO.jl")
include("Mesh/Mesh.jl")
include("P4est/P4est.jl")
include("Parallel/Parallel.jl")
include("Physical_space/Physical_space.jl")
include("Solver/Solver.jl")
include("Theory/Theory.jl")
include("Velocity_space/Velocity_space.jl")

# Runtime initialization: runs once per process when the module is loaded.
# OpenBLAS thread count is process-global runtime state and must be set here
# (not at module top-level, which only runs during precompilation). Each MPI
# rank already provides parallelism, so BLAS stays single-threaded to avoid
# oversubscription / busy-wait spinning.
function __init__()
    BLAS.set_num_threads(1)
end

end # module
