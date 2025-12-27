struct ConfigureForSave{DIM,NDF}
    geometry::Vector{Float64}
    trees_num::Vector{Int64}
    quadrature::Union{Vector{Float64},AbstractQuadrature}
    vs_trees_num::Vector{Int64}
    IC::AbstractInitCondType
    domain::Vector{Domain}
    IB::Vector{AbstractBoundary}
    gas::Gas
    solver::Solver
    user_defined::UDF
end
struct StatusForSave
    max_vs_num::Int
    gradmax::Vector{Float64}
    Δt::Float64
    Δt_ξ::Float64
    sim_time::Float64
    ps_adapt_step::Int
    vs_adapt_step::Int
    partition_step::Int
end
function ConfigureForSave(config::Configure{DIM,NDF}) where{DIM,NDF}
    return ConfigureForSave{DIM,NDF}(
        config.geometry,config.trees_num,config.quadrature,
        config.vs_trees_num,config.IC,config.domain,
        config.IB,config.gas,
        config.solver,
        config.user_defined
    )
end
function StatusForSave(status::Status)
    return StatusForSave(status.max_vs_num,status.gradmax,status.Δt,
        status.Δt_ξ,status.sim_time,status.ps_adapt_step,
        status.vs_adapt_step,status.partition_step
    )
end
function Status(status::StatusForSave)
    return Status(status.max_vs_num,status.gradmax,status.Δt,
        status.Δt_ξ,status.sim_time,status.ps_adapt_step,
        status.vs_adapt_step,status.partition_step,Residual(DIM),Ref(false))
end
struct SolverSet
    config::ConfigureForSave
    mpi_size::Int
end
struct PS_Solution
    prim::Vector{Float64}
    qf::Vector{Float64}
end
struct Boundary_PS_Solution
    prim::Vector{Float64}
    qf::Vector{Float64}
    p::Vector{Float64}
end
function PS_Solution(ps_data::PS_Data{DIM}) where{DIM}
    if ps_data.bound_enc<0
        prim = Vector{Float64}(undef,DIM+2);prim.=NaN
        qf = Vector{Float64}(undef,DIM);qf.=NaN
        return PS_Solution(prim, qf)
    end
    return PS_Solution(ps_data.prim, ps_data.qf)
end
function PS_Solution(::InsideSolidData{DIM,NDF}) where{DIM,NDF}
    prim = Vector{Float64}(undef,DIM+2);prim.=NaN
    qf = Vector{Float64}(undef,DIM);qf.=NaN
    return PS_Solution(prim, qf)
end
struct VS_Solution
    quadid::Int
    ps_midpoint::Vector{Float64}
    midpoint::Matrix{Float64}
    level::Vector{Int8}
    df::Matrix{Float64}
end
struct Solution
    ps_solutions::Vector{PS_Solution}
    vs_solutions::Vector{VS_Solution}
end
function Solution(ps_solutions::Vector{PS_Solution})
    return Solution(ps_solutions, VS_Solution[])
end
struct MeshInfo
    neighbor_nums::Vector{Vector{Int}}
end
struct Result
    solution::Solution
    mesh_info::MeshInfo
end
struct Boundary_Solution
    midpoints::Vector{Vector{Float64}}
    normal::Vector{Vector{Float64}}
    ps_solutions::Vector{Boundary_PS_Solution}
end