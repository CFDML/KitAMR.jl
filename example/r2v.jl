using KitAMR
KitAMR.result2vtk("result2025-05-28_22-41","test_slope_vg")
KitAMR.result2vtk("result2025-05-29_19-26","test_1e-4_level4")
KitAMR.result2vtk("result2025-01-11_22-34","test_vs_vtk")
KitAMR.result2vtk("result2025-01-10_21-28","test_non-vs_vtk")
KitAMR.result2vtk("result2025-01-11_21-29","test_non-adapt_vtk")
KitAMR.result2vtk("result2025-01-14_16-09","test_DVM_vtk")
KitAMR.result2vtk("v100","v100")
KitAMR.result2vtk("result2025-03-06_14-33","test upwind solid corner")
KitAMR.result2vtk("result2025-03-31_17-28","v28")#11600steps with 1TOLERANCE
KitAMR.result2vtk("result2025-02-21_20-27","test adaptive temp")#7000steps with 100TOLERANCE
KitAMR.result2vtk("result2025-05-06_21-43","naca0012")
KitAMR.result2vtk("result2025-04-19_16-14","test_cylinder_vsrefine")
KitAMR.boundary_result2csv("result2025-05-29_19-26","test_1e-4_level4")   