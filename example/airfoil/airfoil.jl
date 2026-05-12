using KitAMR,MPI
include("./airfoil_udf.jl")
MPI.Init()

solver = Solver(;
    DIM = 2, NDF = 2,
    AMR_PS_MAXLEVEL = 7,
    AMR_VS_MAXLEVEL = 0,
    PS_DYNAMIC_AMR = false,
    VS_DYNAMIC_AMR = false,
    flux = CAIDVM,
    time_marching = CAIDVM_Marching,
)
gas = Gas(;
    K = 1.0,
    Kn = 0.026,
    ω = 0.81,
    ωᵣ = 0.81,
)
output = Output(
    solver;
)
udf = UDF(;
    static_ps_refine_flag = shock_wave_region
)
config = Configure(solver;
    geometry = [-3.,7.,-7.,7.],
    trees_num = [24,32],
    quadrature = [-4.,8.,-6.,6.],
    vs_trees_num = [60,60],
    IC = Uniform([1.,2.0*sqrt(5/6),0.,1.]),
    domain = [
            Domain(SuperSonicInflow,1,[1.,2.0*sqrt(5/6),0.,1.]),
            Domain(InterpolatedOutflow,2),Domain(InterpolatedOutflow,3),Domain(InterpolatedOutflow,4)
        ],
    IB = [Vertices(Maxwellian,"./example/airfoil/naca0012.csv",true,4.5,[1.,0.,0.,161/290])],
    output = output,
    gas = gas,
    user_defined = udf
)

p4est,ka = initialize(config)
listen_for_save!()
max_sim_time = 20.
nt = max_sim_time/ka.kinfo.status.Δt+1.0 |> floor |> Int
for i in 1:nt
    adaptive_mesh_refinement!(p4est,ka;partition_interval=160)
    update_slope!(ka)
    slope_exchange!(p4est, ka) 
    update_solid_cell!(ka)
    solid_exchange!(p4est, ka)
    update_solid_neighbor!(ka)
    flux!(ka) 
    iterate!(ka) 
    data_exchange!(p4est, ka)
    check_for_convergence(ka)&&break
    check!(p4est,ka)
end
save_result(p4est,ka)
finalize!(p4est,ka)
MPI.Finalize()

