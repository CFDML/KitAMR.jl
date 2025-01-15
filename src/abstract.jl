const COMM_NUMS_TAG = 10000
const COMM_DATA_TAG = 20000
const EPS = 1e-12
const ADAPT_COEFFI_PS = 0.02
const ADAPT_COEFFI_VS = 1e-4
const SOLID_CELL_ID_NUM = 8
const BALANCE_FLAG = 127
pxest_t = Union{p4est_t,p8est_t}
pxest_ghost_t = Union{p4est_ghost_t,p8est_ghost_t}
pxest_mesh_t = Union{p4est_mesh_t,p8est_mesh_t}

P_pxest_t = Union{Ptr{p4est_t}, Ptr{p8est_t}}
P_pxest_ghost_t = Union{Ptr{p4est_ghost_t}, Ptr{p8est_ghost_t}}
P_pxest_mesh_t = Union{Ptr{p4est_mesh_t}, Ptr{p8est_mesh_t}}
P_pxest_quadrant_t = Union{Ptr{p4est_quadrant_t}, Ptr{p8est_quadrant_t}}
P_pxest_iter_volume_info_t = Union{Ptr{p4est_iter_volume_info_t}, Ptr{p8est_iter_volume_info_t}}
P_pxest_iter_face_info_t = Union{Ptr{p4est_iter_face_info_t}, Ptr{p8est_iter_face_info_t}}

PW_pxest_t = Union{PointerWrapper{p4est_t}, PointerWrapper{p8est_t}}
PW_pxest_ghost_t = Union{PointerWrapper{p4est_ghost_t}, PointerWrapper{p8est_ghost_t}}
PW_pxest_mesh_t = Union{PointerWrapper{p4est_mesh_t}, PointerWrapper{p8est_mesh_t}}
PW_pxest_quadrant_t = Union{PointerWrapper{p4est_quadrant_t}, PointerWrapper{p8est_quadrant_t}}
PW_pxest_iter_volume_info_t = Union{PointerWrapper{p4est_iter_volume_info_t}, PointerWrapper{p8est_iter_volume_info_t}}
PW_pxest_iter_face_info_t = Union{PointerWrapper{p4est_iter_face_info_t},PointerWrapper{p8est_iter_face_info_t}}