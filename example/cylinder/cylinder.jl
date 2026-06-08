using KitAMR,MPI
include("./cylinder_udf.jl")
MPI.Init()

solver = Solver(;
    DIM = 2, NDF = 2,
    AMR_PS_MAXLEVEL = 7,
    AMR_DYNAMIC_PS_MAXLEVEL = 4,
    AMR_VS_MAXLEVEL = 3,
    PS_DYNAMIC_AMR = true,
    VS_DYNAMIC_AMR = true,
    flux = CAIDVM,
    time_marching = CAIDVM_Marching,
    max_sim_time = 20.,
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
    # static_ps_refine_flag = shock_wave_region
    dynamic_ps_refine_flag = amr_region
)
config = Configure(solver;
    geometry = [-16.,16.,-16.,16.],
    trees_num = [25,25],
    quadrature = [-10.,10.,-10.,10.],
    vs_trees_num = [16,16],
    IC = PCoordFn(cylinder_buffer_IC),
    domain = [
            Domain(SuperSonicInflow,1,[1.,5.0*sqrt(5/6),0.,1.]),Domain(UniformOutflow,2),
            Domain(UniformOutflow,3),Domain(UniformOutflow,4)
        ],
    IB = [Circle(Maxwellian,[0.,0.],1.,true,4.0,[1.,0.,0.,1.])],
    output = output,
    gas = gas,
    user_defined = udf
)


p4est,ka = initialize(config);
solve!(p4est, ka; prerefine_steps = 0, ps_interval = 40, partition_interval = 40)
save_result(p4est,ka)
finalize!(p4est,ka)
MPI.Finalize()

