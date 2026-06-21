include("common.jl")

RP3DScaling.with_mpi() do
    RP3DScaling.run_strong_initial_amr(64)
end
