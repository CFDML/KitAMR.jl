using KitAMR,MPI
MPI.Init()

solver = Solver(;
    DIM = 3, NDF = 1,
    AMR_PS_MAXLEVEL = 0,
    AMR_VS_MAXLEVEL = 0,
    PS_DYNAMIC_AMR = false,
    VS_DYNAMIC_AMR = false,
    flux = CAIDVM,
    time_marching = CAIDVM_Marching,
)
gas = Gas(;
    K = 0.0,
    Kn = 0.075,
    ω = 0.81,
    ωᵣ = 0.81,
)
output = Output(solver)
udf = UDF()
config = Configure(solver;
    geometry = [-0.5,0.5,-0.5,0.5,-0.5,0.5],
    trees_num = [8,8,64],
    quadrature = [-5.0,5.0,-5.0,5.0,-5.0,5.0],
    vs_trees_num = [24,24,24],
    IC = Uniform([1.,0.,0.,0.,1.]),
    domain = [
            Domain(Period,1),
            Domain(Period,2),Domain(Period,3),Domain(Period,4),
            Domain(Maxwellian,5,[1.,0.,0.,0.,1.]),
            Domain(Maxwellian,6,[1.,1.0*sqrt(5/6),0.,0.,1.0])
        ],
    output = output,
    gas = gas,
    user_defined = udf
)

p4est,ka = initialize(config);
max_sim_time = 20.
nt = max_sim_time/ka.kinfo.status.Δt+1.0 |> floor |> Int
for i in 1:10
    # adaptive_mesh_refinement!(p4est,ka)
    if MPI.Comm_rank(MPI.COMM_WORLD)==0
        @show i     
        if i == 5
            time_start = time_ns()  
        end
        if i==10
            time_end = time_ns()
            global time_start
            time = (time_end-time_start)/1e9
            @show time
        end
    end
    update_slope!(ka)
    slope_exchange!(p4est, ka) 
    update_solid_cell!(ka)
    solid_exchange!(p4est, ka)
    update_solid_neighbor!(ka)
    flux!(ka) 
    iterate!(ka) 
    data_exchange!(p4est, ka)
    check_for_convergence(ka)&&break
end
# save_result(p4est,ka)
finalize!(p4est,ka)
MPI.Finalize()

