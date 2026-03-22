using KitAMR, MPI
MPI.Init()

solver = Solver(;
    DIM = 2, NDF = 2,
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
    geometry = [-0.5,0.5],
    trees_num = [16,16],
    quadrature = [-5.0,5.0,-5.0,5.0],
    vs_trees_num = [16,16],
    IC = Uniform([1.,0.,0.,1.]),
    domain = [
            Domain(Maxwellian,1,[1.,0.,0.,1.]),
            Domain(Maxwellian,2,[1.,1.0*sqrt(5/6),0.,1.0]),Domain(Period,3),
            Domain(Period,4)
        ],
    output = output,
    gas = gas,
    user_defined = udf
)

ps4est,amr = initialize_KitAMR(config);
max_sim_time = 20.
nt = max_sim_time/amr.global_data.status.Δt+1.0 |> floor |> Int
for i in 1:10
    # adaptive_mesh_refinement!(ps4est,amr)
    update_slope!(amr)
    slope_exchange!(ps4est, amr) 
    update_solid_cell!(amr)
    solid_exchange!(ps4est, amr)
    update_solid_neighbor!(amr)
    flux!(amr) 
    iterate!(amr) 
    data_exchange!(ps4est, amr)
    check_for_convergence(amr)&&break
    check!(i,ps4est,amr)
end
# save_result(ps4est,amr)
finalize!(ps4est,amr)
MPI.Finalize()
