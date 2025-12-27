using KitAMR
include("convergence_ps_TC/convergence_ps_udf.jl")
include("artificial_interpolation/convergence_ps_udf.jl")
include("cylinder/cylinder_udf.jl")
KitAMR.result2vtk("result2025-01-09_21-43","modified1_UGKS_vtk")
KitAMR.result2vtk("result2025-08-17_19-27","test_period")
KitAMR.result2vtk("result2025-01-11_22-34","test_vs_vtk")
KitAMR.result2vtk("result2025-01-10_21-28","test_non-vs_vtk")
KitAMR.result2vtk("result2025-01-11_21-29","test_non-adapt_vtk")
KitAMR.result2vtk("result2025-01-14_16-09","test_DVM_vtk")
KitAMR.result2vtk("v100","v100")
KitAMR.result2vtk("result2025-07-25_18-16","hypersonic_cylinder_vuniform_vanleer")
KitAMR.result2vtk("result2025-08-24_16-59","test cylinder_Kn0p01 downwind_reduce")#11600steps with 1TOLERANCE
KitAMR.result2vtk("result2025-06-17_11-43","test cylinder_conserved")#7000steps with 100TOLERANCE
KitAMR.result2vtk("result2025-07-11_06-03","naca0012_vuniform")
KitAMR.result2vtk("result2025-07-10_23-44","hypersonic_cylinder_vuniform")
for i in [64,128,256]
    dir = "result_interp_p"*string(i)
    name = "interp_p"*string(i)
    KitAMR.boundary_result2csv(dir,name)       
end
KitAMR.boundary_result2csv("result2025-12-08_23-36","cylinder_vanleer")   