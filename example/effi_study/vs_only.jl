include("common.jl")

R2DVJEffi.with_mpi() do
    R2DVJEffi.run_mode("vs_only")
end
