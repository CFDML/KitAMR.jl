const COMM_NUMS_TAG = 10000
const COMM_DATA_TAG = 20000
const EPS = 1e-12
const ADAPT_COEFFI_PS = 0.25
const ADAPT_COEFFI_VS = 0.125
const SOLID_CELL_ID_NUM = 8
const BALANCE_FLAG = 127
const RES_CHECK_INTERVAL = 100
const TOLERANCE = 1e-6
const REDUNDANT_STEPS_NUM = 100
const TIME_STEP_CONTRACT_RATIO = 0.2
const face_num_2d = 4
const face_num_3d = 6
const WLST = [6,10]
const RFT = [[],[[2, 4], [1, 4], [3, 2], [1, 3]],[[-1,-1,-1],[1,-1,-1],[-1,1,-1],[1,1,-1],[-1,-1,1],[1,-1,1],[-1,1,1],[1,1,1]]] # refine face table    
const RMT = [[],[[-1, -1], [1, -1], [-1, 1], [1, 1]],[[-1,-1,-1],[1,-1,-1],[-1,1,-1],[1,1,-1],[-1,-1,1],[1,-1,1],[-1,1,1],[1,1,1]]] # refine midpoint table
const ANTIVT = [[],[[-1, -1], [1, -1], [1, 1], [-1, 1]],[[-1,-1,-1],[1,-1,-1],[-1,1,-1],[1,1,-1],[-1,-1,1],[1,-1,1],[-1,1,1],[1,1,1]]] # anticlock vertices table
const NMT = [[],[[-1,0],[1,0],[0,-1],[0,1]],[[-1,0,0],[1,0,0],[0,-1,0],[0,1,0],[0,0,-1],[0,0,1]]] # neighbor midpoint table
const FAT = [[2,1],[[2,3],[1,3],[1,2]]]# face area table
const CLP = [1,3,5,7] # Vertices' indices in cut_rect
const CTN = [[],[[1,3],[2,3],[1,4],[2,4]],[]] # Corner target neighbor DIM{Corner_id{possible faces}}
const pxest_ts = [p4est_t,p8est_t]
const pxest_quadrant_ts = [p4est_quadrant_t,p8est_quadrant_t]
const pxest_ghost_ts = [p4est_ghost_t,p8est_ghost_t]
const pxest_mesh_ts = [p4est_mesh_t,p8est_mesh_t]
const pxest_iter_volume_info_ts = [p4est_iter_volume_info_t,p8est_iter_volume_info_t]
const pxest_iter_face_info_ts = [p4est_iter_face_info_t,p8est_iter_face_info_t]
const pxest_iter_face_side_ts = [p4est_iter_face_side_t,p8est_iter_face_side_t]

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
PW_pxest_iter_face_side_t = Union{PointerWrapper{p4est_iter_face_side_t},PointerWrapper{p8est_iter_face_side_t}}

"""
$(TYPEDEF)
"""
abstract type AbstractInitCondType end

abstract type AbstractVTKCellType end
"""
$(TYPEDEF)
"""
abstract type Pixel <: AbstractVTKCellType end
"""
$(TYPEDEF)
"""
abstract type Triangle <: AbstractVTKCellType end
"""
$(TYPEDEF)
"""
abstract type Voxel <: AbstractVTKCellType end
"""
$(TYPEDEF)
"""
abstract type Tetra <: AbstractVTKCellType end

"""
$(TYPEDEF)
"""
abstract type AbstractFluxType end
"""
$(TYPEDEF)
"""
abstract type AbstractDVMFluxType <: AbstractFluxType end
"""
$(TYPEDEF)
"""
abstract type AbstractTimeMarchingType end


abstract type UGKS<:AbstractDVMFluxType end
abstract type DVM<:AbstractDVMFluxType end
"""
$(TYPEDEF)
Conserved adaptive implicit DVM flux.
"""
abstract type CAIDVM<:AbstractDVMFluxType end # Conserved DVM, which encures energy conserved flux.
const MicroFlux = Union{DVM}
const HybridFlux = Union{UGKS,CAIDVM} # to add more...


abstract type Euler <:AbstractTimeMarchingType end
abstract type UGKS_Marching<:AbstractTimeMarchingType end
"""
$(TYPEDEF)
"""
abstract type CAIDVM_Marching<:AbstractTimeMarchingType end
"""
$(TYPEDEF)
"""
abstract type AbstractPsData{DIM,NDF} end
abstract type AbstractGhostPsData{DIM,NDF} <: AbstractPsData{DIM,NDF} end

"""
$(TYPEDEF)
"""
abstract type AbstractFace end
abstract type InnerFace <: AbstractFace end
abstract type BoundaryFace <: AbstractFace end

abstract type AbstractGas end

"""
$(TYPEDEF)
"""
abstract type AbstractBoundaryType end
"""
$(TYPEDEF)
"""
abstract type AbstractBoundCondType end

"""
$(TYPEDEF)
Fully diffused Maxwellian gas-surface interaction model.
"""
abstract type Maxwellian<:AbstractBoundCondType end

"""
$(TYPEDEF)
"""
abstract type SuperSonicInflow <: AbstractBoundCondType end
"""
$(TYPEDEF)
"""
abstract type SuperSonicOutflow <: AbstractBoundCondType end
"""
$(TYPEDEF)
"""
abstract type UniformOutflow <: AbstractBoundCondType end
"""
$(TYPEDEF)
"""
abstract type InterpolatedOutflow <: AbstractBoundCondType end
abstract type AxisSymmetric <: AbstractBoundCondType end
"""
$(TYPEDEF)
"""
abstract type Period <: AbstractBoundCondType end
const AbstractBCType = Union{Vector,Function}

"""
$(TYPEDEF)
"""
abstract type AbstractVsData{DIM,NDF} end

abstract type AbstractQuadrature end

include("../Gas/Types.jl")
include("../Velocity_space/Types.jl")
include("../Physical_space/Types.jl")
include("../Boundary/Types.jl")
include("../Solver/Types.jl")
include("../IO/Types.jl")

export SuperSonicInflow, SuperSonicOutflow, UniformOutflow, InterpolatedOutflow, Period, Maxwellian
export AbstractInitCondType, AbstractBoundCondType, AbstractBoundaryType
export Pixel, Voxel, Triangle, Tetra
export AbstractPsData, AbstractFace
export AbstractVsData
export AbstractFluxType, AbstractDVMFluxType, CAIDVM