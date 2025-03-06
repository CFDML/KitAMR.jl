using KitAMR
KitAMR.result2vtk("result2025-01-09_21-43","modified1_UGKS_vtk")
KitAMR.result2vtk("result2025-01-14_15-20","test_UGKS_vtk")
KitAMR.result2vtk("result2025-01-11_22-34","test_vs_vtk")
KitAMR.result2vtk("result2025-01-10_21-28","test_non-vs_vtk")
KitAMR.result2vtk("result2025-01-11_21-29","test_non-adapt_vtk")
KitAMR.result2vtk("result2025-01-14_16-09","test_DVM_vtk")
KitAMR.result2vtk("result2025-02-24_16-02","test_newIB")
KitAMR.result2vtk("result2025-03-06_14-33","test upwind solid corner")
KitAMR.result2vtk("result2025-02-21_22-32","test vanleer IB")#11600steps with 1TOLERANCE
KitAMR.result2vtk("result2025-02-21_20-27","test adaptive temp")#7000steps with 100TOLERANCE
KitAMR.result2vtk("result2025-02-24_20-05","./MixedBoundary/grid")
KitAMR.result2vtk("result2025-02-20_08-41","./CylinderMa5Kn0p1Ar/static refine/grid")
KitAMR.boundary_result2csv("result2025-03-06_14-33","./test upwind solid step back")