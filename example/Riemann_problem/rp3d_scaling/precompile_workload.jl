include("common.jl")

RP3DScaling.with_mpi() do
    RP3DScaling.run_precompile_workload()
end
