#using Pkg
#Pkg.activate(@__DIR__)
using MPI
using P4est
#using P4estTypes
include("../lib/P4estTypes/src/P4estTypes.jl")
using .P4estTypes

using StaticArrays
using Parameters
using SpecialFunctions
#using Makie
using CairoMakie
using PythonCall
scipy = pyimport("scipy")
np = pyimport("numpy")
using LinearAlgebra
using HDF5
using JLD2
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
# include("flux.jl")
include("flux_interp_t.jl")
include("iterate.jl")
include("partition.jl")
include("postprocess.jl")
include("finalize.jl")
