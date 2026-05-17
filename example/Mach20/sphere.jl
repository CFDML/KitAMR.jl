using KitAMR,MPI
include("./sphere_udf.jl")
MPI.Init()

solver = Solver(;
    DIM = 3, NDF = 1,
    CFL = 0.2,
    AMR_PS_MAXLEVEL = 4,
    AMR_DYNAMIC_PS_MAXLEVEL = 3,
    AMR_VS_MAXLEVEL = 5,
    PS_DYNAMIC_AMR = true,
    VS_DYNAMIC_AMR = true,
    flux = CAIDVM,
    time_marching = CIP_Marching,
    ADAPT_COEFFI_PS = 0.5,
    ADAPT_COEFFI_VS_LOCAL = 0.05,
    ADAPT_COEFFI_VS_GLOBAL = 0.25
)
gas = Gas(;
    K = 0.0,
    Kn = 1.0,
    ω = 0.75,
    ωᵣ = 0.75,
    μᵣ = 5.0*√π/16.0*1.0
)
output = Output(solver)
udf = UDF(;
    # dynamic_ps_refine_flag = shock_wave_region
)
config = Configure(solver;
    geometry = [-4.,4.,-4.,4.,-4.,4.],
    trees_num = [16,16,16],
    quadrature = [-57.94,57.94,-57.94,57.94,-57.94,57.94], # 3σ: 5√Ts
    vs_trees_num = [8,8,8],
    IC = PCoordFn(sphere_buffer_IC),
    domain = [
        Domain(SuperSonicInflow,1,[1.,20.0*sqrt(5/6),0.,0.,1.]),
        Domain(UniformOutflow,2),Domain(UniformOutflow,3),Domain(UniformOutflow,4),
        Domain(UniformOutflow,5),Domain(UniformOutflow,6)
    ],
    IB = [Sphere(Maxwellian,[0.,0.,0.],0.5,true,1.5,[1.,0.,0.,0.,1/(1.0+(5/3-1)*0.5*20.0^2)])],
    output = output,
    gas = gas,
    user_defined = udf
)

p4est,ka = initialize(config);
KitAMR.execute_check(p4est,ka)
listen_for_save!()
max_sim_time = 20.
nt = max_sim_time/ka.kinfo.status.Δt+1.0 |> floor |> Int
for i in 1:nt
    if MPI.Comm_rank(MPI.COMM_WORLD)==0
        @show i
    end
    adaptive_mesh_refinement!(p4est,ka;ps_interval = 20, vs_interval = 20, partition_interval = 40, vs_balance = false)
    slope!(p4est,ka) # Update `sdf` in `AbstractVsData` and `sw` in `PsData`.
    flux!(p4est, ka) # Compute and update numerical fluxes.
    iterate!(p4est, ka) # Collision process and time marching.
    check!(p4est,ka) # Check for save and output simulation status to `stdout`.
    check_for_convergence(ka)&&break # Check for convergence.
end
KitAMR.execute_check(p4est,ka)
save_result(p4est,ka)
finalize!(p4est,ka)
MPI.Finalize()

