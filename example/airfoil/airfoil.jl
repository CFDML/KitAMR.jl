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
    max_sim_time = 20.,
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

p4est,ka = initialize(config; prerefine_steps = 0)
solve!(p4est, ka; partition_interval = 160)
save_result(p4est,ka)
finalize!(p4est,ka)
MPI.Finalize()

