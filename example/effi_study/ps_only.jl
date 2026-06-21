include("common.jl")

R2DVJEffi.with_mpi() do
    R2DVJEffi.run_mode("ps_only")
end
