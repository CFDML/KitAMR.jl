module P4estTypes

using AbstractTrees: print_tree
using DocStringExtensions
using MPI
using P4est
using RecipesBase
using SparseArrays

# export QuadrantWrapper, Connectivity, Pxest, GhostLayer, LNodes

# export refine, brick, pxest
# export lengthoflocalquadrants, lengthofglobalquadrants # Find better names
# export level, unsafe_storeuserdata!, unsafe_loaduserdata
# export offset, coordinates
# export iterateforest
# export refine!, coarsen!, balance!, partition!, expand!
# export lnodes, ghostlayer, connectivity
# export ghosts, mirrors, sharers
# export globalid
export Connectivity

# Running
#     grep -C 2 "closures are not supported" src/codegen.cpp in the julia
# source will show what archs do not support cfunction closures.  The julia
# tracking issue from AArch64 is
#    https://github.com/JuliaLang/julia/issues/27174
const cfunction_closure = Sys.ARCH === :x86_64 || Sys.ARCH === :i686

include("sc.jl")

initialized() = P4est.package_id() >= 0
function initialize()
    if !SC.initialized()
        sc_init(MPI.COMM_NULL, 0, 0, C_NULL, SC_LP_ERROR)
    end
    if !initialized()
        p4est_init(C_NULL, SC_LP_ERROR)
    end
end

"""
    P4estTypes.setverbosity(logpriority::P4estTypes.SC.LP.LogPriority)

Sets the verbosity of p4est. See [`P4estTypes.SC.LP.LogPriority`](@ref P4estTypes.SC.LP.LogPriority)
for a list of valid priorities.
"""
function setverbosity(logpriority::SC.LP.LogPriority)
    sc_package_set_verbosity(P4est.package_id(), logpriority)
end

"""
    P4estTypes.uses_mpi()

Returns `true` if the p4est C library uses MPI.
"""
uses_mpi() = P4est.uses_mpi()

include("connectivity.jl")
include("pxest.jl")
include("lnodes.jl")
include("ghost.jl")

function __init__()
    initialize()
end

using PrecompileTools
@setup_workload begin
    @compile_workload begin
        MPI.Initialized() || MPI.Init()
        initialize()
        p = pxest(brick(3, 4))
        p = pxest(brick(3, 4, 2))
    end
end

end # module
