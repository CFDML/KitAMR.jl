using KitAMR,MPI
include("./cylinder_udf.jl")
MPI.Init()

solver = Solver(;
    DIM = 2, NDF = 2,
    AMR_PS_MAXLEVEL = 7,
    AMR_PS_DYNAMIC_MAXLEVEL = 5,
    AMR_VS_MAXLEVEL = 4,
    AMR_PS_DYNAMIC = true,
    AMR_VS_DYNAMIC = true,
    AMR_VS_MODE = :haar,
    AMR_PS_THRES = 0.2,
    flux = CAIDVM,
    time_marching = CIP_Marching,
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
    vtk_celltype = [Pixel,Triangle],
)
udf = UDF(;
    dynamic_ps_refine_flag = amr_region
)
config = Configure(solver;
    geometry = [-16.,16.,-16.,16.],
    trees_num = [16,16],
    quadrature = [-19.839,19.839,-19.839,19.839],
    vs_trees_num = [8,8],
    IC = PCoordFn(cylinder_buffer_IC),
    domain = [
            Domain(SuperSonicInflow,1,[1.,5.0*sqrt(5/6),0.,1.]),Domain(UniformOutflow,2),
            Domain(UniformOutflow,3),Domain(UniformOutflow,4)
        ],
    IB = [Circle(Maxwellian,[0.,0.],1.,true,4.,[1.,0.,0.,1.])],
    output = output,
    gas = gas,
    user_defined = udf
)


p4est,ka = initialize(config);
solve!(p4est, ka;vs_balance=true,ps_interval = 40, vs_interval = 40, partition_interval = 40)
save_result(p4est,ka)
save_for_restart(p4est,ka;dir_path="cylinder_restart")
finalize!(p4est,ka)
MPI.Finalize()

