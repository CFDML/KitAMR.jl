include("AMR.jl")
include("Criteria.jl")
include("Cut_cell.jl")
include("Initialize.jl")
export VS_Data, Ghost_VS_Data, Face_VS_Data, CuttedVelocityCells
export vs_adaptive_mesh_refinement!, vs_refine!, vs_coarsen!, vs_conserved_correction!, contribution_refine_flag, contribution_coarsen_flag