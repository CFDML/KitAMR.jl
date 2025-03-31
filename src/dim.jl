const face_num_2d = 4
const face_num_3d = 6

const WLST = [6,10]
const RFT = [[],[[2, 4], [1, 4], [3, 2], [1, 3]],[[-1,-1,-1],[1,-1,-1],[-1,1,-1],[1,1,-1],[-1,-1,1],[1,-1,1],[-1,1,1],[1,1,1]]] # refine face table    
const RMT = [[],[[-1, -1], [1, -1], [-1, 1], [1, 1]],[[-1,-1,-1],[1,-1,-1],[-1,1,-1],[1,1,-1],[-1,-1,1],[1,-1,1],[-1,1,1],[1,1,1]]] # refine midpoint table
const ANTIVT = [[],[[-1, -1], [1, -1], [1, 1], [-1, 1]],[[-1,-1,-1],[1,-1,-1],[-1,1,-1],[1,1,-1],[-1,-1,1],[1,-1,1],[-1,1,1],[1,1,1]]] # anticlock vertices table
const NMT = [[],[[-1,0],[1,0],[0,-1],[0,1]]] # neighbor midpoint table
const FAT = [[2,1],[[2,3],[1,3],[1,2]]]# face area table
const CLP = [1,3,5,7] # Vertices' indices in cut_rect
const pxest_ts = [p4est_t,p8est_t]
const pxest_quadrant_ts = [p4est_quadrant_t,p8est_quadrant_t]
const pxest_ghost_ts = [p4est_ghost_t,p8est_ghost_t]
const pxest_mesh_ts = [p4est_mesh_t,p8est_mesh_t]
const pxest_iter_volume_info_ts = [p4est_iter_volume_info_t,p8est_iter_volume_info_t]
const pxest_iter_face_info_ts = [p4est_iter_face_info_t,p8est_iter_face_info_t]
const pxest_iter_face_side_ts = [p4est_iter_face_side_t,p8est_iter_face_side_t]




# if DIM===2
#     #=
#     2D Z-ordering
#     y
#     |-----|-----|
#     |  3  |  4  |
#     |-----|-----|
#     |  1  |  2  |
#     |-----|-----| x
#     =#
#     const moment_u = moment_u_2D
#     const micro_slope = micro_slope_2D
#     const rft = [[2, 4], [1, 4], [3, 2], [1, 3]] # refine face table    
#     const rmt = [[-1, -1], [1, -1], [-1, 1], [1, 1]] # refine midpoint table
#     const cft = [1, 2, 1, 3]
#     const nft = [[1, 3], [2, 4], [1, 2], [3, 4]] # face neighbor inner-sided table, faceid{z-ordering}
#     const AMR_CONNECT_FACE = P4EST_CONNECT_FACE
#     const AMR_ghost_new = p4est_ghost_new
#     const AMR_mesh_new_ext = p4est_mesh_new_ext
# end
# if DIM === 3
#     #=
#     3D Z-ordering
#             +--------+--------+
#             /        /        /|
#     z      /    7   /    8   / |
#     |     /        /        /  |
#         +--------+--------+  8|
#         /        /        / | /|
#     /    5   /    6   /  |/ |
#     /        /        / 6 /  |
#     +--------+--------+   /| 4|
#     |        |        |  / | /
#     |   5    |    6   | /  |/ 
#     |        |        |/   /  
#     +--------+--------+ 2 /   
#     |        |        |  /    
#     |    1   |    2   | /     
#     |        |        |/      
#     +--------+--------+   -x
#     =#
#     # const moment_u = moment_u_3D
#     # const micro_slope = micro_slope_3D
#     const rmt3 = [[-1,-1,-1],[1,-1,-1],[-1,1,-1],[1,1,-1],[-1,-1,1],[1,-1,1],[-1,1,1],[1,1,1]]
#     const AMR_CONNECT_FACE = P8EST_CONNECT_FACE
#     const AMR_ghost_new = p8est_ghost_new
#     const AMR_mesh_new_ext = p8est_mesh_new_ext
# end
# if DIM === 2 && NDF === 2
#     const discrete_maxwell = discrete_maxwell_2DF_2D
#     const shakhov_part = shakhov_part_2DF_2D
#     const calc_fwb = calc_fwb_2DF_2D
#     const calc_w0 = calc_w0_2DF_2D
#     const calc_flux_f0 = calc_flux_f0_2DF_2D
#     const calc_qf = calc_qf_2DF_2D
# end
# if DIM === 3 && NDF === 2
#     # const discrete_maxwell = discrete_maxwell_2DF_3D
#     # const shakhov_part = shakhov_part_2DF_3D
#     # const calc_fwb = calc_fwb_2DF_3D
#     # const calc_w0 = calc_w0_2DF_3D
#     # const calc_flux_f0 = calc_flux_f0_2DF_3D
#     # const calc_qf = calc_qf_2DF_3D
# end