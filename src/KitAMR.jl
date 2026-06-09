module KitAMR

using MPI
using JLD2
using Dates
using LinearAlgebra
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
using ProgressMeter

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
    # Quiet p4est's C library (libsc): by default it prints an informational line for every
    # forest operation (new/refine/coarsen/balance/partition/ghost). Raise its log threshold
    # so only genuine errors are shown. Override with the `KITAMR_P4EST_LOG` environment
    # variable, set to an `SC_LP_*` integer threshold — e.g. `6` (SC_LP_PRODUCTION) to restore
    # the per-operation chatter, `8` (SC_LP_ERROR, the default here), `9` (SC_LP_SILENT) to mute
    # everything including errors.
    threshold = something(tryparse(Cint, get(ENV, "KITAMR_P4EST_LOG", "")), Cint(SC_LP_ERROR))
    p4est_init(C_NULL, threshold)
end

end # module
