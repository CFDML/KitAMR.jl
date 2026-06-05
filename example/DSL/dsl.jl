using KitAMR,MPI
MPI.Init()
function dsl_init(midpoint,kinfo)
    Ma = 2
    u0 = √(5/6)*Ma
    k = 80; δ = 0.05
    T0 = 5/6
    uy = u0*δ*sin(2*π*(midpoint[1]+1/4))
    if midpoint[2]<0.5
        ux = u0*tanh(k*(midpoint[2]-1/4))
        return [1.,ux,uy,1/T0]
    else
        ux = u0*tanh(k*(3/4-midpoint[2]))
        return [1.,ux,uy,1/T0]
    end
end
solver = Solver(;
    DIM = 2, NDF = 2,
    CFL = 0.4,
    AMR_PS_MAXLEVEL = 3,
    AMR_VS_MAXLEVEL = 3,
    PS_DYNAMIC_AMR = true,
    ADAPT_COEFFI_PS = 0.3,
    VS_DYNAMIC_AMR = true,
    flux = CAIDVM,
    # flux = DVM,
    time_marching = CIP_Marching,
    # time_marching = Euler,
    # time_marching = CAIDVM_Marching,
    max_sim_time = 20.,
)
gas = Gas(;
    K = 1.0,
    Kn = 1e-5,
    ω = 0.81,
    ωᵣ = 0.81,
)
output = Output(
    solver;
)
udf = UDF(;
)
config = Configure(solver;
    geometry = [0., 1., 0., 1.],
    trees_num = [16,16],
    quadrature = [-7.,7.,-7.,7.],
    vs_trees_num = [8,8],
    IC = PCoordFn(dsl_init),
    domain = [
            Domain(Period,1),Domain(Period,2),
            Domain(Period,3),Domain(Period,4)
        ],
    output = output,
    gas = gas,
    user_defined = udf
)

p4est,ka = initialize(config);
slope!(p4est,ka)
for _ in 1:3
    KitAMR.ps_adaptive_mesh_refinement!(p4est,ka;recursive = false)
    KitAMR.amr_recover!(p4est,ka)
end
KitAMR.execute_check!(p4est,ka)
listen_for_save!()
check_for_animsave!(p4est, ka)   # initial animation frame (no-op unless output.anim_dt>0)
i = 0
while !reached_max_time(ka)
    global i += 1
    if MPI.Comm_rank(MPI.COMM_WORLD)==0
        @show i
    end
    adaptive_mesh_refinement!(p4est,ka;ps_interval = 20, vs_interval = 40, partition_interval = 20, vs_balance = false)
    limit_Δt!(ka) # Shrink Δt to land exactly on the next animation frame / max_sim_time.
    slope!(p4est,ka) # Update `sdf` in `AbstractVsData` and `sw` in `PsData`.
    flux!(p4est, ka) # Compute and update numerical fluxes.
    iterate!(p4est, ka) # Collision process and time marching.
    check_for_animsave!(p4est, ka) # Write a frame if this step landed on a frame time.
    check!(p4est,ka) # Check for save and output simulation status to `stdout`.
    check_for_convergence(ka)&&break # Check for convergence.
end
KitAMR.execute_check!(p4est,ka)
save_result(p4est,ka)
finalize!(p4est,ka)
MPI.Finalize()
