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
solve!(p4est, ka;
    prerefine_steps = 3,
    ps_interval = 20, vs_interval = 40, partition_interval = 20)
save_result(p4est,ka)
finalize!(p4est,ka)
MPI.Finalize()
