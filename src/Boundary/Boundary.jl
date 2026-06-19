include("AMR.jl")
include("Circle.jl")
include("Immersed_boundary.jl")
include("Parallel.jl")
include("Vertices.jl")
include("Period.jl")
include("Positivity.jl")
include("Triangles.jl")

export Domain, CompositeBC, Circle, Sphere, Vertices, Triangles, TriangleKDT
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

# Composite weights are always converted to a convex blend. This keeps cases such as
# `(a, 1-a)` and `(100a, 100(1-a))` equivalent and rejects ambiguous negative weights.
function _normalize_composite_weights(raw_weights, n::Integer)
    weights = Float64.(collect(raw_weights))
    length(weights) == n ||
        error("CompositeBC returned $(length(weights)) weights for $n component domains.")
    @inbounds for weight in weights
        isfinite(weight) || error("CompositeBC weights must be finite; got $weight.")
        weight >= 0.0 || error("CompositeBC weights must be non-negative; got $weight.")
    end
    total = sum(weights)
    total > 0.0 || error("CompositeBC weights must contain at least one positive entry.")
    weights ./= total
    return weights
end

# Runtime flux evaluation can use the full face/solver context when a weight function needs it.
function composite_weights(weights::AbstractVector, n::Integer, face::DomainFace, ka::KA)
    return _normalize_composite_weights(weights, n)
end
function composite_weights(weights::Function, n::Integer, face::DomainFace, ka::KA)
    midpoint = face.midpoint
    raw = if applicable(weights, face, ka)
        weights(face, ka)
    elseif applicable(weights, midpoint, ka.kinfo)
        weights(midpoint, ka.kinfo)
    elseif applicable(weights, midpoint)
        weights(midpoint)
    else
        weights(; midpoint)
    end
    return _normalize_composite_weights(raw, n)
end
composite_weights(bc::CompositeBC, face::DomainFace, ka::KA) =
    composite_weights(bc.weights, length(bc.domains), face, ka)

# AMR and initial velocity-space refinement only have a boundary midpoint and `KInfo`.
function composite_weights(weights::AbstractVector, n::Integer,
                           midpoint::AbstractVector, kinfo::KInfo)
    return _normalize_composite_weights(weights, n)
end
function composite_weights(weights::Function, n::Integer,
                           midpoint::AbstractVector, kinfo::KInfo)
    raw = if applicable(weights, midpoint, kinfo)
        weights(midpoint, kinfo)
    elseif applicable(weights, midpoint)
        weights(midpoint)
    else
        weights(; midpoint)
    end
    return _normalize_composite_weights(raw, n)
end
composite_weights(bc::CompositeBC, midpoint::AbstractVector, kinfo::KInfo) =
    composite_weights(bc.weights, length(bc.domains), midpoint, kinfo)

# Boundary primitive sampling is shared by Lohner/auto-AMR sensors and initialization. Composite
# components without a primitive state, for example `UniformOutflow`, are ignored here because
# their effect is already represented by the interior state during flux evaluation.
function domain_bc_prim(domain::Domain, midpoint::AbstractVector, kinfo::KInfo)
    isdefined(domain, :bc) || return nothing
    return domain_bc_prim(domain.bc, midpoint, kinfo)
end
domain_bc_prim(bc::AbstractVector, midpoint::AbstractVector, kinfo::KInfo) =
    Float64.(collect(bc))
function domain_bc_prim(bc::Function, midpoint::AbstractVector, kinfo::KInfo)
    if applicable(bc, midpoint, kinfo)
        return Float64.(collect(bc(midpoint, kinfo)))
    elseif applicable(bc, midpoint)
        return Float64.(collect(bc(midpoint)))
    end
    return Float64.(collect(bc(; midpoint)))
end
function domain_bc_prim(bc::CompositeBC, midpoint::AbstractVector, kinfo::KInfo)
    weights = composite_weights(bc, midpoint, kinfo)
    prim = nothing
    prim_weight = 0.0
    @inbounds for i in eachindex(bc.domains)
        weights[i] == 0.0 && continue
        component_prim = domain_bc_prim(bc.domains[i], midpoint, kinfo)
        component_prim === nothing && continue
        prim_weight += weights[i]
        if prim === nothing
            prim = weights[i] .* component_prim
        else
            length(component_prim) == length(prim) ||
                error("CompositeBC component primitive length mismatch.")
            prim .+= weights[i] .* component_prim
        end
    end
    prim === nothing && return nothing
    prim ./= prim_weight
    return prim
end
