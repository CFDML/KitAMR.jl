using KitAMR,MPI
MPI.Init()
function Sod_init(midpoint,kinfo)
    if midpoint[1]<0.
        return [1.,10√(5/6),0.,0.1]
    else
        return [10.,0.,0.,1.]
    end
end
solver = Solver(;
    DIM = 2, NDF = 2,
    CFL = 0.2,
    AMR_PS_MAXLEVEL = 0,
    AMR_DYNAMIC_PS_MAXLEVEL = 0,
    AMR_VS_MAXLEVEL = 5,
    PS_DYNAMIC_AMR = false,
    VS_DYNAMIC_AMR = true,
    flux = CAIDVM,
    # flux = DVM,
    time_marching = CIP_Marching,
    # time_marching = Euler,
    # time_marching = CAIDVM_Marching,
)
gas = Gas(;
    K = 1.0,
    Kn = 0.1,
    ω = 0.81,
    ωᵣ = 0.81,
)
output = Output(
    solver;
)
udf = UDF(;
)
config = Configure(solver;
    geometry = [-2.,2.,-0.5,0.5],
    trees_num = [128,32],
    quadrature = [-41.43,41.43,-41.43,41.43],
    vs_trees_num = [8,8],
    IC = PCoordFn(Sod_init),
    domain = [
            Domain(UniformOutflow,1),Domain(UniformOutflow,2),
            Domain(Period,3),Domain(Period,4)
        ],
    output = output,
    gas = gas,
    user_defined = udf
)

p4est,ka = initialize(config);
slope!(p4est,ka)
KitAMR.ps_adaptive_mesh_refinement!(p4est,ka;recursive = true)
KitAMR.execute_check(p4est,ka)
listen_for_save!()
max_sim_time = 20.
nt = max_sim_time/ka.kinfo.status.Δt+1.0 |> floor |> Int
for i in 1:5000
    adaptive_mesh_refinement!(p4est,ka;ps_interval = 20, vs_interval = 10, partition_interval = 20, vs_balance = true)
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
