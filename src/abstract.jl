const AV = AbstractVector
const AM = AbstractMatrix
const AA = AbstractArray
const AVOM = AbstractVecOrMat
const IN = Integer
const RN = Real
const FN = AbstractFloat
abstract type AbstractControlVolume end
abstract type AbstractInterface end
abstract type AbstractVelocitySpace end
abstract type AbstractPhysicalSpace end
abstract type AbstractGas end
abstract type AbstractCondition end
abstract type AbstractInitial end
abstract type AbstractBoundary end
abstract type AbstractSolution end
abstract type AbstractPhysicalSpace2D <: AbstractPhysicalSpace end
abstract type AbstractVelocitySpace2D <: AbstractVelocitySpace end
abstract type AbstractFlux end
abstract type AbstractFlux2F <: AbstractFlux end
const DIM = 2
const NDF = 2
# const GHOST_OFFSET = DIM + NDF
const DVM_VS_MAXLEVEL = 3
const DVM_PS_MAXLEVEL = 3
const COMM_NUMS_TAG = 10000
const COMM_DATA_TAG = 20000
const EPS = 1e-12
