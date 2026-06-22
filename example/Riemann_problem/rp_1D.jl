using KitAMR,MPI
MPI.Init()
function Sod_init(midpoint,kinfo)
    if midpoint[1]<0.
        return [1.,0.,0.,1.]
    else
        return [0.125,0.,0.,0.625]
    end
end


solver = Solver(;
    DIM = 2, NDF = 2,
    CFL = 0.4,
    AMR_PS_MAXLEVEL = 4,
    AMR_VS_MAXLEVEL = 4,
    AMR_PS_DYNAMIC = true,
    AMR_PS_THRES = 0.25,
    AMR_VS_DYNAMIC = true,
    flux = CAIDVM,
    time_marching = CIP_Marching,
    max_sim_time = 0.5,
)
gas = Gas(;
    K = 1.0,
    Kn = 1e3,
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
solve!(p4est, ka)
save_result(p4est,ka)
finalize!(p4est,ka)
MPI.Finalize()
