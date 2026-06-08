using KitAMR,MPI
include("./X38_udf.jl")

MPI.Init() # MPI initialization. Mandatory for a paralleled program using MPI.



# ----------------------------------------------------------------------------------------------------
#=
Configuration by directly construct `Configure` struct.
=#
solver = Solver(;
    DIM = 3, NDF = 1,
    AMR_PS_MAXLEVEL = 4,
    AMR_DYNAMIC_PS_MAXLEVEL = 4,
    AMR_VS_MAXLEVEL = 3,
    PS_DYNAMIC_AMR = true,
    VS_DYNAMIC_AMR = true,
    flux = CAIDVM,
    time_marching = CAIDVM_Marching,
    max_sim_time = 20.,
)
gas = Gas(;
    K = 0.0,
    Kn = 0.275,
    ω = 0.81,
    ωᵣ = 0.81,
    T_ref = 273/56,
)
output = Output(solver)
udf = UDF(;
    dynamic_ps_refine_flag = amr_region
)
config = Configure(solver;
    geometry = [-4.0,5.0,-4.0,4.0,-4.0,4.0],
    trees_num = [16,16,16],
    quadrature = [-14.2,21.3,-17.75,17.75,-17.75,17.75],
    vs_trees_num = [10,10,10],
    IC = PCoordFn(X38_buffer_IC),
    domain = [
            Domain(SuperSonicInflow,1,
            [1.,8.0*sqrt(5/6),0.,0.,1.0]),Domain(UniformOutflow,2),Domain(UniformOutflow,3),
            Domain(UniformOutflow,4),Domain(UniformOutflow,5),Domain(UniformOutflow,6)
        ],
    IB = [Triangles(Maxwellian,"./example/X38/X38_normalized.stl",true,1.5,[1.,0.,0.,0.,56/300])],
    output = output,
    gas = gas,
    user_defined = udf
)
# ----------------------------------------------------------------------------------------------------

p4est,ka = initialize(config) # Initialization for `KitAMR_Data`.
#=
Run the time-stepping loop. `max_sim_time` is set in the `Solver` above; `solve!`
terminates (and lands the last step) exactly on it via `reached_max_time` / `limit_Δt!`.
=#
solve!(p4est, ka;
    prerefine_steps = 0,
    ps_interval = 40, vs_interval = 40, partition_interval = 40)
save_result(p4est,ka) # Save converging results.
finalize!(p4est,ka) # Finalize `p4est` things. Release the memory managed by `C`.
MPI.Finalize() # MPI finalization. Mandatory for a paralleled program using MPI.

