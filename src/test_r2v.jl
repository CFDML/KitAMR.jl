using KitAMR
KitAMR.result2vtk("result2025-01-09_21-43","modified1_UGKS_vtk")
KitAMR.result2vtk("result2025-01-14_15-20","test_UGKS_vtk")
KitAMR.result2vtk("result2025-01-11_22-34","test_vs_vtk")
KitAMR.result2vtk("result2025-01-10_21-28","test_non-vs_vtk")
KitAMR.result2vtk("result2025-01-11_21-29","test_non-adapt_vtk")
KitAMR.result2vtk("result2025-01-14_16-09","test_DVM_vtk")
KitAMR.result2vtk("result2025-01-15_15-57","test_abs_vtk")
KitAMR.result2vtk("result2025-02-14_16-30","static_test")
KitAMR.result2vtk("result2025-02-15_22-46","testR05")#11600steps with 1TOLERANCE
KitAMR.result2vtk("result2025-02-10_15-09","Mach5_Kn001")#7000steps with 100TOLERANCE
KitAMR.result2vtk("result2025-02-15_22-12","./CylinderMach5Kn01/CylinderMach5Kn01_filtered")
KitAMR.result2vtk("result2025-02-14_18-47","Mach5_Kn01_Static")
KitAMR.boundary_result2csv("result2025-02-15_22-12","./CylinderMach5Kn01/boundary_filtered")