struct SolverSet
    config::Configure
    mpi_size::Int
end
struct PS_Solution
    prim::Vector{Float64}
    qf::Vector{Float64}
end
function PS_Solution(ps_data::PS_Data)
    if ps_data.bound_enc<0
        prim = Vector{Float64}(undef,DIM+2);prim.=NaN
        qf = Vector{Float64}(undef,DIM);qf.=NaN
        return PS_Solution(prim, qf)
    end
    return PS_Solution(ps_data.prim, ps_data.qf)
end
function PS_Solution(ps_data::InsideSolidData{DIM,NDF}) where{DIM,NDF}
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