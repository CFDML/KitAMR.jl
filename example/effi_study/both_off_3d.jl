ENV["VJE_DIM"] = "3"
include("common.jl")

R2DVJEffi.with_mpi() do
    R2DVJEffi.run_mode("both_off")
end
