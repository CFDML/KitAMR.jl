struct ConfigureForSave{DIM,NDF}
    geometry::Vector{Float64}
    trees_num::Vector{Int64}
    quadrature::Vector{Float64}
    vs_trees_num::Vector{Int64}
    IC::AbstractICType
    domain::Vector{Domain}
    IB::Vector{AbstractBoundary}
    gas::Gas
    solver::Solver
end
function ConfigureForSave(config::Configure{DIM,NDF}) where{DIM,NDF}
    return ConfigureForSave{DIM,NDF}(
        config.geometry,config.trees_num,config.quadrature,
        config.vs_trees_num,config.IC,config.domain,
        config.IB,config.gas,
        config.solver
    )
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