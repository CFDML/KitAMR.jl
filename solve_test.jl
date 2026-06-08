# Validates the new solve! user API on two configs:
#   (A) no immersed boundary  (B) a circle immersed boundary.
# Each: initialize + solve! capped by max_steps; checks it runs and returns ka.
# Run: mpirun -np 2 julia --project=. solve_test.jl
using KitAMR, MPI
MPI.Init()
RANK = MPI.Comm_rank(MPI.COMM_WORLD)

function run_case(name, config)
    p4est, ka = initialize(config)
    # finalize=false so we can inspect ka afterwards; save/animation off for a quiet test.
    solve!(p4est, ka;
        ps_interval = 5, vs_interval = 10, partition_interval = 5,
        max_steps = 20,
        listen_for_save = false, animation = false, status_check = false,
        )
    step = ka.kinfo.status.step
    RANK == 0 && println("[$name] solve! ran to step = $step")
    KitAMR.finalize!(p4est, ka)
    return step
end

# ---- (A) no immersed boundary ----
solverA = Solver(; DIM=2, NDF=2, AMR_PS_MAXLEVEL=2, AMR_DYNAMIC_PS_MAXLEVEL=2, AMR_VS_MAXLEVEL=1,
                   PS_DYNAMIC_AMR=true, VS_DYNAMIC_AMR=true, flux=CAIDVM, time_marching=CAIDVM_Marching)
configA = Configure(solverA;
    geometry=[-0.5,0.5,-0.5,0.5], trees_num=[8,8], quadrature=[-5.0,5.0,-5.0,5.0], vs_trees_num=[8,8],
    IC=Uniform([1.,0.,0.,1.]),
    domain=[Domain(Maxwellian,1,[1.,0.,0.,1.]),Domain(Maxwellian,2,[1.,1.0*sqrt(5/6),0.,1.0]),
            Domain(Period,3),Domain(Period,4)],
    output=Output(solverA), gas=Gas(; K=0.0,Kn=0.075,ω=0.81,ωᵣ=0.81), user_defined=UDF())
sA = run_case("no-IB", configA)

# ---- (B) circle immersed boundary ----
solverB = Solver(; DIM=2, NDF=2, AMR_PS_MAXLEVEL=3, AMR_DYNAMIC_PS_MAXLEVEL=2, AMR_VS_MAXLEVEL=1,
                   PS_DYNAMIC_AMR=true, VS_DYNAMIC_AMR=true, flux=CAIDVM, time_marching=CAIDVM_Marching)
configB = Configure(solverB;
    geometry=[-8.,8.,-8.,8.], trees_num=[8,8], quadrature=[-5.0,5.0,-5.0,5.0], vs_trees_num=[8,8],
    IC=Uniform([1.,0.,0.,1.]),
    domain=[Domain(SuperSonicInflow,1,[1.,2.0,0.,1.]),Domain(UniformOutflow,2),
            Domain(UniformOutflow,3),Domain(UniformOutflow,4)],
    IB=[Circle(Maxwellian,[0.,0.],1.,true,4.0,[1.,0.,0.,1.])],
    output=Output(solverB), gas=Gas(; K=0.0,Kn=0.075,ω=0.81,ωᵣ=0.81), user_defined=UDF())
sB = run_case("circle-IB", configB)

if RANK == 0
    ok = sA > 0 && sB > 0
    println(ok ? "SOLVE TEST OK" : "SOLVE TEST FAILED")
end
MPI.Finalize()
