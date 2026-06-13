using KitAMR,MPI
include("./sphere_udf.jl")
MPI.Init()

solver = Solver(;
    DIM = 3, NDF = 1,
    AMR_PS_MAXLEVEL = 4,
    AMR_DYNAMIC_PS_MAXLEVEL = 3,
    AMR_VS_MAXLEVEL = 3,
    PS_DYNAMIC_AMR = true,
    VS_DYNAMIC_AMR = true,
    flux = CAIDVM,
    time_marching = CAIDVM_Marching,
    ADAPT_COEFFI_PS = 0.5,
    ADAPT_COEFFI_VS_LOCAL = 0.05,
    ADAPT_COEFFI_VS_GLOBAL = 0.25,
    max_sim_time = 20.,
)
gas = Gas(;
    K = 0.0,
    Kn = 0.1,
    ω = 0.75,
    ωᵣ = 0.75,
    μᵣ = 5.0*√π/16.0*0.1
)
output = Output(solver)
udf = UDF(;
    dynamic_ps_refine_flag = shock_wave_region
)
config = Configure(solver;
    geometry = [-4.,4.,-4.,4.,-4.,4.],
    trees_num = [16,16,16],
    quadrature = [-12.14,12.14,-12.14,12.14,-12.14,12.14], # 3σ: 3√Ts
    vs_trees_num = [8,8,8],
    IC = PCoordFn(sphere_buffer_IC),
    domain = [
        Domain(SuperSonicInflow,1,[1.,3.834*sqrt(5/6),0.,0.,1.]),
        Domain(UniformOutflow,2),Domain(UniformOutflow,3),Domain(UniformOutflow,4),
        Domain(UniformOutflow,5),Domain(UniformOutflow,6)
    ],
    IB = [Sphere(Maxwellian,[0.,0.,0.],0.5,true,1.5,[1.,0.,0.,0., 1.0/(1.0+(5/3-1)*0.5*3.834^2)])],
    output = output,
    gas = gas,
    user_defined = udf
)

p4est,ka = initialize(config; prerefine_steps = 0);
solve!(p4est, ka;
    ps_interval = 40, vs_interval = 40, partition_interval = 40)
save_result(p4est,ka)
finalize!(p4est,ka)
MPI.Finalize()

