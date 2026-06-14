include("AMR.jl")
include("Circle.jl")
include("Immersed_boundary.jl")
include("Parallel.jl")
include("Vertices.jl")
include("Period.jl")
include("Positivity.jl")
include("Triangles.jl")

export Domain, Circle, Sphere, Vertices, Triangles, TriangleKDT
export DomainFace, ImmersedBoundary
export initialize_solid_neighbor!, update_solid_neighbor!, update_solid_cell!
export initialize_immersed_boundaries!
export solid_exchange!, solid_exchange_begin!, solid_exchange_finish!


function overlap_test(lower,upper,hyper_rec::HyperRectangle{SVector{N,Float64}}) where N
    for i = 1:N
        (upper[i]<hyper_rec.mins[i]||lower[i]>hyper_rec.maxes[i])&& return false
    end
    return true
end

get_bc(bc::AbstractVector; kwargs...) = collect(bc)
function get_bc(bc::Function; kwargs...)
    if haskey(kwargs, :intersect_point) && haskey(kwargs, :ib)
        intersect_point = kwargs[:intersect_point]
        ib = kwargs[:ib]
        applicable(bc, intersect_point, ib) && return collect(bc(intersect_point, ib))
        return collect(bc(; intersect_point, ib))
    elseif haskey(kwargs, :midpoint)
        midpoint = kwargs[:midpoint]
        applicable(bc, midpoint) && return collect(bc(midpoint))
        return collect(bc(; midpoint))
    end
    error("Function boundary condition requires `midpoint` for domain boundaries or `intersect_point` and `ib` for immersed boundaries.")
end

function get_domain_bc(face::DomainFace, ka::KA)
    get_domain_bc(face.domain.bc, face, ka)
end
function get_domain_bc(bc::AbstractVector, face::DomainFace, ka::KA)
    get_bc(bc)
end
function get_domain_bc(bc::Function, face::DomainFace, ka::KA)
    midpoint = face.midpoint
    if applicable(bc, midpoint, ka.kinfo)
        return collect(bc(midpoint, ka.kinfo))
    elseif applicable(bc, midpoint)
        return collect(bc(midpoint))
    end
    collect(bc(; midpoint))
end
