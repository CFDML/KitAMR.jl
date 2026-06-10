# Domain boundary
"""
$(TYPEDEF)

Structure of domain boundary.

## Fields

$(TYPEDFIELDS)
"""
struct Domain{T<:AbstractBoundCond} <: AbstractBoundary
    """
    Index of the domain boundary. From 1 to 6, it represents the boundary of `xmin`, `xmax`, `ymin`, `ymax`, `zmin`, `zmax`.
    """
    id::Int
    """
    Whether refine at the domain boundary. Default is `false`.
    """
    refine::Bool
    """
    The boundary condition at the domain boundary.
    """
    bc::AbstractBCType
    Domain(T, id; refine = false) = new{T}(id, refine)
    Domain(T, id, bc; refine = false) = new{T}(id, refine, bc)
end

"""
$(TYPEDEF)
Face at the domain edge.
"""
struct DomainFace{DIM,NDF,T}<:BoundaryFace
    rot::Float64
    direction::Int
    midpoint::Vector{Float64}
    domain::Domain{T}
    ps_data::PsData{DIM,NDF}
end

# Immersed boundary
"""
$(TYPEDEF)
$(TYPEDFIELDS)
"""
struct Circle{T<:AbstractBoundCond} <: AbstractBoundary
    """
    Center of the circle.
    """
    center::Vector
    """
    Radius of the circle.
    """
    radius::Real
    """
    Is solid inside the circle? # If is solid inside the circle, it should be true. Otherwise, it should be false.
    """
    solid::Bool # If is solid inside the circle, it should be true. Otherwise, it should be false.
    """
    Refinement coefficient. # The ratio of the search radius and the minimal mesh scale, which determines the maximal distance between IB nodes and aux_points.
    """
    search_coeffi::Real # The ratio of the search radius and the minimal mesh scale, which determines the maximal distance between IB nodes and aux_points.
    """
    Primary macroscopic variables of the solid.
    """
    bc::AbstractBCType
    """
    Maximum distance of the refinement region from the boundary.
    """
    search_radius::Real
    Circle(
        ::Type{T},
        center::Vector,
        radius,
        solid,
        search_coeffi,
        bc,
    ) where {T<:AbstractBoundCond} = new{T}(center, radius, solid, search_coeffi, bc)
    Circle(c::Circle{T}, ds::Float64) where {T<:AbstractBoundCond} =
        new{T}(c.center, c.radius, c.solid, c.search_coeffi, c.bc, c.search_coeffi*ds)
end
"""
$(TYPEDEF)
$(TYPEDFIELDS)
"""
struct Sphere{T<:AbstractBoundCond} <: AbstractBoundary # 3D circle
    center::Vector
    radius::Real
    solid::Bool
    search_coeffi::Real
    bc::AbstractBCType
    search_radius::Real
    Sphere(
        ::Type{T},
        center::Vector,
        radius,
        solid,
        search_coeffi,
        bc,
    ) where {T<:AbstractBoundCond} = new{T}(center, radius, solid, search_coeffi, bc)
    Sphere(c::Sphere{T}, ds::Float64) where {T<:AbstractBoundCond} =
        new{T}(c.center, c.radius, c.solid, c.search_coeffi, c.bc, c.search_coeffi*ds)
end
const AbstractCircle = Union{Circle,Sphere}

"""
$(TYPEDEF)
$(TYPEDFIELDS)
"""
struct Vertices{DIM,T<:AbstractBoundCond} <: AbstractBoundary
    """
    Vertices of the boundary, sorted in clockwise or counterclockwise order.
    """
    vertices::Vector{Vector{Float64}}
    solid::Bool # Is solid inside the boundary?
    bc::AbstractBCType
    """
    The outer box of the vertices as `[[xmin,ymin,zmin],[xmax,ymax,zmax]]`.
    """
    box::Vector{Vector{Float64}}
    """
    The number of cell layers refined from the boundary.
    """
    search_radius::Real
end
"""
$(TYPEDSIGNATURES)

  - `file` is the path to the `.csv` file.
"""
function Vertices(
    ::Type{T},
    file::String,
    solid,
    refine_coeffi,
    bc,
) where {T<:AbstractBoundCond}
    s = CSV.read(file, DataFrame; header = true)
    DIM = length(names(s))
    vertices = [[s.x[i], s.y[i]] for i in eachindex(s.x)]
    if vertices[1]==vertices[end]
        vertices = vertices[1:(end-1)]
    end
    box = [[minimum(s.x), minimum(s.y)], [maximum(s.x), maximum(s.y)]]
    return Vertices{DIM,T}(vertices, solid, bc, box, refine_coeffi)
end
function Vertices(v::Vertices{DIM,T}, config) where {DIM,T<:AbstractBoundCond}
    ds_max = maximum([
        (config[:geometry][2i]-config[:geometry][2i-1])/config[:trees_num][i] for
        i = 1:config[:DIM]
    ])
    ds = norm([
        (config[:geometry][2i]-config[:geometry][2i-1])/config[:trees_num][i]/2^config[:AMR_PS_MAXLEVEL]
        for i = 1:config[:DIM]
    ])
    return Vertices{DIM,T}(
        v.vertices,
        v.solid,
        v.bc,
        [v.box[1] .- ds_max, v.box[2] .+ ds_max],
        v.search_radius*ds,
    )
end

"""
$(TYPEDEF)
$(TYPEDFIELDS)
"""
struct TriangleKDT
    kdt::KDTree
    mesh::Mesh
    table::Vector{HyperRectangle{SVector{3,Float64}}}
    triangle_recs::Vector{HyperRectangle{SVector{3,Float64}}} # Accessed by id of mesh.faces
    triangle_edges::Vector{Vector{SVector{3,Float64}}} # Accessed by id of mesh.faces
end
function TriangleKDT(mesh::Mesh)
    centers = zeros(3, length(mesh.faces))
    for i in eachindex(mesh.faces)
        vertices = mesh.position[mesh.faces[i]]
        for j = 1:3
            for k = 1:3
                centers[j, i] += vertices[k][j]/3.0
            end
        end
    end
    kdt = KDTree(centers; leafsize = 24)
    return TriangleKDT(
        kdt,
        mesh,
        triangle_box_table(kdt, mesh),
        triangle_recs(mesh),
        triangle_edges(mesh),
    )
end

"""
$(TYPEDEF)
$(TYPEDFIELDS)
"""
struct Triangles{T<:AbstractBoundCond} <: AbstractBoundary
    solid::Bool
    bc::AbstractBCType
    search_radius::Real
    tkdt::TriangleKDT
end
"""
$(TYPEDSIGNATURES)

  - `file` is the path to the `.stl` file.
"""
function Triangles(
    ::Type{T},
    file::String,
    solid,
    search_radius,
    bc,
) where {T<:AbstractBoundCond}
    mesh = load(file)
    tkdt = TriangleKDT(mesh)
    return Triangles{T}(solid, bc, search_radius, tkdt)
end
function Triangles(
    ::Type{T},
    solid::Bool,
    search_radius,
    bc,
    tkdt::TriangleKDT,
    config,
) where {T<:AbstractBoundCond}
    ds = norm([
        (config[:geometry][2i]-config[:geometry][2i-1])/config[:trees_num][i]/2^config[:AMR_PS_MAXLEVEL]
        for i = 1:config[:DIM]
    ])
    return Triangles{T}(solid, bc, search_radius*ds, tkdt)
end
function Triangles(ib::Triangles{T}, config) where {T}
    return Triangles(T, ib.solid, ib.search_radius, ib.bc, ib.tkdt, config)
end



mutable struct SolidCells{DIM,NDF}
    ps_datas::Vector{PsData{DIM,NDF}}
    quadids::Vector{Cint} # A better choice: store all quadid of SolidCells to avoid extensive comparing iteration. Are only needed to be updated before partition.
    #=
    quadids provide the information of solidcells on other processors. Only after this, the exchange of IB nodes can be executeable.
    =#
end
function SolidCells(ps_datas::Vector{PsData{DIM,NDF}}) where {DIM,NDF}
    quadids = [x.quadid for x in ps_datas]
    return SolidCells{DIM,NDF}(ps_datas, quadids)
end

"""
$(TYPEDEF)
Pre-collected immersed-boundary data for a simulation.  Donor cells, solid
cells and IB-adjacent faces are gathered once at initialisation (and after
AMR) so that the time-step loop can iterate over them directly without
scanning every cell in the tree.
$(TYPEDFIELDS)
"""
mutable struct ImmersedBoundary{DIM,NDF}
    """
    Donor cells (`bound_enc > 0`).
    """
    donor_cells::Vector{AbstractPsData{DIM,NDF}}
    """
    Solid cells (`bound_enc < 0`, excluding `InsideSolidData`).
    """
    solid_cells::Vector{AbstractPsData{DIM,NDF}}
    """
    Faces that involve a `SolidNeighbor` on one side.
    """
    faces::Vector{AbstractFace}
end
