using KitAMR
KitAMR.result2vtk("result2025-01-09_21-43","modified1_UGKS_vtk")
KitAMR.result2vtk("result2025-01-14_15-20","test_UGKS_vtk")
KitAMR.result2vtk("result2025-01-11_22-34","test_vs_vtk")
KitAMR.result2vtk("result2025-01-10_21-28","test_non-vs_vtk")
KitAMR.result2vtk("result2025-01-11_21-29","test_non-adapt_vtk")
KitAMR.result2vtk("result2025-01-14_16-09","test_DVM_vtk")
KitAMR.result2vtk("v100","v100")
KitAMR.result2vtk("result2025-03-06_14-33","test upwind solid corner")
KitAMR.result2vtk("result2025-06-24_20-11","test cylinder_Kn0p01")#11600steps with 1TOLERANCE
KitAMR.result2vtk("result2025-06-17_11-43","test cylinder_conserved")#7000steps with 100TOLERANCE
KitAMR.result2vtk("result2025-06-04_22-17","naca0012_conserved")
KitAMR.result2vtk("result2025-06-21_16-35","test_cylinder")
KitAMR.boundary_result2csv("result2025-06-21_10-05","cylinder_Kn0p01")   