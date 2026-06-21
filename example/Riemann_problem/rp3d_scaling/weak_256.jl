include("common.jl")

RP3DScaling.with_mpi() do
    RP3DScaling.run_weak_uniform(256)
end
