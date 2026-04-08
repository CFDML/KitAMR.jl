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

p4est,ka = initialize_KitAMR(config) # Initialization for `KitAMR_Data`.
listen_for_save!() # Start listening for `save` input from `stdin`.
max_sim_time = 20. # Maximum simulation time.
nt = max_sim_time/ka.kinfo.status.Δt+1.0 |> floor |> Int # Maximum number of time marching steps.


#=
Main loop.
=#
for i in 1:nt
    adaptive_mesh_refinement!(p4est,ka;ps_interval = 40, vs_interval=40, partition_interval=40) # AMR.
    update_slope!(ka) # Update `sdf` in `VsData` and `sw` in `PsData`.
    slope_exchange!(p4est, ka) # Update `sdf` in `GhostVsData` by MPI communication.
    update_solid_cell!(ka) # Update variables in solid cells in immersed boundaries.
    solid_exchange!(p4est, ka) # Update variables in ghost solid cells by MPI communication.
    update_solid_neighbor!(ka) # Update variables in SolidNeighbor by IBM.
    flux!(ka) # Compute numerical fluxes.
    iterate!(ka) # Collision process and time marching.
    data_exchange!(p4est, ka) # Update variables in ghost cells by MPI communication.
    check_for_convergence(ka)&&break # Check for convergence.
    check!(p4est,ka) # Check for save and output simulation status to `stdout`.
end

save_result(p4est,ka) # Save converging results.
finalize!(p4est,ka) # Finalize `p4est` things. Release the memory managed by `C`.
MPI.Finalize() # MPI finalization. Mandatory for a paralleled program using MPI.

