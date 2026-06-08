using KitAMR,MPI
MPI.Init()
function Sod_init(midpoint,kinfo)
    if midpoint[1]<0.
        return [1.,0.,0.,0.5]
    else
        return [0.125,0.,0.,0.625]
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
    # time_marching = CIP_Marching,
    time_marching = Euler,
    # time_marching = CAIDVM_Marching,
    max_sim_time = 20.,
)
gas = Gas(;
    K = 1.0,
    Kn = 0.001,
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
    trees_num = [32,8],
    quadrature = [-7.,7.,-7.,7.],
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
solve!(p4est, ka;
    prerefine_steps = 1, prerefine_recursive = true,
    ps_interval = 20, vs_interval = 20, partition_interval = 20)
save_result(p4est,ka)
finalize!(p4est,ka)
MPI.Finalize()
