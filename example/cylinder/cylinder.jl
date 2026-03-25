using KitAMR,MPI
include("./cylinder_udf.jl")
MPI.Init()

solver = Solver(;
    DIM = 2, NDF = 2,
    AMR_PS_MAXLEVEL = 7,
    AMR_VS_MAXLEVEL = 0,
    PS_DYNAMIC_AMR = false,
    VS_DYNAMIC_AMR = false,
    flux = CAIDVM,
    time_marching = CAIDVM_Marching,
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
    static_ps_refine_flag = shock_wave_region
)
config = Configure(solver;
    geometry = [-16.,16.,-16.,16.],
    trees_num = [25,25],
    quadrature = [-10.,10.,-10.,10.],
    vs_trees_num = [60,60],
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


p4est,ka = initialize_KitAMR(config);
listen_for_save!()
max_sim_time = 20.
nt = max_sim_time/ka.kinfo.status.Δt+1.0 |> floor |> Int
for i in 1:nt
    adaptive_mesh_refinement!(p4est,ka;partition_interval=160)
    update_slope!(ka)
    slope_exchange!(p4est, ka) 
    update_solid_cell!(ka)
    solid_exchange!(p4est, ka)
    update_solid_neighbor!(ka)
    slope_exchange!(p4est, ka) 
    flux!(ka) 
    iterate!(ka) 
    data_exchange!(p4est, ka)
    check_for_convergence(ka)&&break
    check!(p4est,ka)
    # KitAMR.check_for_animsave!(p4est,ka)
    # if ka.kinfo.status.sim_time>3.0
    #     break
    # end
end
save_result(p4est,ka)
finalize!(p4est,ka)
MPI.Finalize()

