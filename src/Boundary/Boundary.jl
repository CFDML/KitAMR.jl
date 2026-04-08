include("AMR.jl")
include("Circle.jl")
include("Immersed_boundary.jl")
include("Parallel.jl")
include("Vertices.jl")
include("Period.jl")
include("Triangles.jl")
get_bc(bc::AbstractVector;kwargs...) = copy(bc)
function get_bc(bc::Function;intersect_point,ib,kwargs...)
    bc(;intersect_point,ib)
end

export Domain, Circle, Sphere, Vertices, Triangles, TriangleKDT
export DomainFace
export initialize_solid_neighbor!, update_solid_neighbor!, update_solid_cell!
export solid_exchange!