include("AMR.jl")
include("Circle.jl")
include("Immersed_boundary.jl")
include("Parallel.jl")
include("Vertices.jl")
include("Period.jl")
include("Triangles.jl")

export Domain, Circle, Sphere, Vertices, Triangles, TriangleKDT
export DomainFace
export initialize_solid_neighbor!, update_solid_neighbor!, update_solid_cell!
export solid_exchange!


function overlap_test(lower,upper,hyper_rec::HyperRectangle{SVector{N,Float64}}) where N
    for i = 1:N
        (upper[i]<hyper_rec.mins[i]||lower[i]>hyper_rec.maxes[i])&& return false
    end
    return true
end

get_bc(bc::AbstractVector;kwargs...) = copy(bc)
function get_bc(bc::Function;intersect_point,ib,kwargs...)
    bc(;intersect_point,ib)
end