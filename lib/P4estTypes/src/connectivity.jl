"""
    Connectivity{X,P}

Connectivity for a `Pxest` which holds the mesh information for
the roots of `Pxest` quadtrees or octrees. The parameter `X` is 4 if the
roots are quads (2D aka p4est) and 8 if they are hexes (3D aka p8est).

# Fields
$(DocStringExtensions.FIELDS)

# Usage

    Connectivity{4}(name::Symbol)

Construct a connectivity mesh for the roots of a forest-of-quadtrees using
p4est's built-in mesh connectivities.  Valid values for `name` are

- `:unitsquare`: the unit square.
- `:periodic`: all-periodic unit square.
- `:rotwrap`: periodic unit square (the left and right faces are identified,
  and bottom and top opposite).
- `:corner`: three-tree mesh around a corner.
- `:pillow`: two trees on top of each other.
- `:moebius`: a five-tree moebius band.
- `:star`: six-tree star.
- `:cubed`: six sides of a unit cube.
- `:disk_nonperiodic`: five-tree flat spherical disk.
- `:icosahedron`: for mapping the sphere using an icosahedron (see
  `@doc P4estTypes.P4est.p4est_connectivity_new_icosahedron` for more info).
- `:shell2d`: 2D spherical shell.
- `:disk2d`: maps a 2D disk.


---
    Connectivity{4}(:disk, periodic_x::Bool, periodic_y::Bool)

Create a connectivity structure for a five-tree flat spherical disk.  The
arguments `periodic_x` and `periodic_y` determine if the disk is periodic in
the x and y directions, respectively.

See `@doc P4estTypes.P4est.p4est_connectivity_new_disk` for detailed information.

---
    Connectivity{8}(name::Symbol)

Construct a connectivity mesh for the roots of a forest-of-octrees using
p8est's built-in mesh connectivities.  Valid values for `name` are

- `:unitcube`: the unit cube.
- `:periodic`: an all-periodic unit cube.
- `:rotcubes`: contains a few cubes (these are rotated against each other to
  stress the topology routines).
- `:rotwrap`: a mostly periodic unit cube (see
  `@doc P4estTypes.P4est.p8est_connectivity_new_rotwrap`).
- `:shell`: a spherical shell (see
  `@doc P4estTypes.P4est.p8est_connectivity_new_shell`).
- `:sphere`: a solid sphere (see
  `@doc P4estTypes.P4est.p8est_connectivity_new_sphere`).
- `:twocubes`: two cubes.
- `:twowrap`: two cubes where the two far ends are identified periodically.

---
    Connectivity{8}(:torus, nsegments)

Create a connectivity structure that builds a revolution torus.  Here
`nsegments` are the number of trees along the great circle.

See `@doc P4estTypes.P4est.p8est_connectivity_new_torus` for detailed
information.

---
    Connectivity{8}(:torus, nsegments)

Create a connectivity structure that builds a revolution torus.  Here
`nsegments` are the number of trees along the great circle.

See `@doc P4estTypes.P4est.p8est_connectivity_new_torus` for detailed
information.

---
    Connectivity{X}(:twotrees, l_face, r_face, orientation) where {X}

Create a connectivity structure (`X=4` for quadtrees and `X=8` for octrees) for
two trees being rotated with respect to each other in a user-defined way.  Here
`l_face` and `r_face` are the 0-based indices of left and right faces,
respectively. The argument `orientation` gives the orientation code of the
trees with respect to each other.

---
    Connectivity{X}(vertices, elements) where {X}

Creates a connectivity from the given list of vertices and element-to-vertex
connectivity.  The parameter set to `X=4` is for quads and `X=8` for hexes.

- `vertices`: should be a number-of-vertices by 3 matrix where the columns
  correspond to x, y, and z coordinates (typically the `z` coordinate will be
  zero for a 2D forest).
- `elements`: should be a number-of-vertices by 4 or 8 matrix where the columns
  vertex indices used to define each element. Note that z-ordering should be
  used, and it should use zero-indexing.

---
    Connectivity{X}(filename::String) where {X}

Create a connectivity from an ABAQUS input at `filename`. The parameter set to
`X=4` is for quads and `X=8` for hexes.

See `@doc P4estTypes.P4est.p4est_connectivity_read_inp` and
`@doc P4estTypes.P4est.p8est_connectivity_read_inp` for example ABAQUS input
files.

# See also
- [`brick`](@ref): a function to create a rectangular [`Connectivity`](@ref).
- [`connectivity`](@ref): a function to get the connectivity of a [`Pxest`](@ref).
- [`refine`](@ref): a function to create a refined [`Connectivity`](@ref).
"""
mutable struct Connectivity{X,P}
    """The pointer (of type `P`) can be a pointer to either a
    `P4estTypes.P4est.p4est_connectivity` or a
    `P4estTypes.P4est.p8est_connectivity`.  See the help
    documentation for these types for more information about the underlying
    p4est structures. """
    pointer::P
    function Connectivity{4}(pointer::Ptr{P4est.LibP4est.p4est_connectivity})
        connectivity = new{4,typeof(pointer)}(pointer)
        # finalizer(connectivity) do c
        #     p4est_connectivity_destroy(c.pointer)
        #     c.pointer = C_NULL
        #     return
        # end
    end
    function Connectivity{8}(pointer::Ptr{P4est.LibP4est.p8est_connectivity})
        connectivity = new{8,typeof(pointer)}(pointer)
        # finalizer(connectivity) do c
        #     p8est_connectivity_destroy(c.pointer)
        #     c.pointer = C_NULL
        #     return
        # end
    end
end

function Connectivity{4}(name::Symbol, args...)
    if name == :unitsquare
        pointer = p4est_connectivity_new_unitsquare()
    elseif name == :periodic
        pointer = p4est_connectivity_new_periodic()
    elseif name == :rotwrap
        pointer = p4est_connectivity_new_rotwrap()
    elseif name == :twotrees
        pointer = p4est_connectivity_new_twotrees(args...)
    elseif name == :corner
        pointer = p4est_connectivity_new_corner()
    elseif name == :pillow
        pointer = p4est_connectivity_new_pillow()
    elseif name == :moebius
        pointer = p4est_connectivity_new_moebius()
    elseif name == :star
        pointer = p4est_connectivity_new_star()
    elseif name == :cubed
        pointer = p4est_connectivity_new_cubed()
    elseif name == :disk_nonperiodic
        pointer = p4est_connectivity_new_disk_nonperiodic()
    elseif name == :disk
        pointer = p4est_connectivity_new_disk(args...)
    elseif name == :icosahedron
        pointer = p4est_connectivity_new_icosahedron()
    elseif name == :shell2d
        pointer = p4est_connectivity_new_shell2d()
    elseif name == :disk2d
        pointer = p4est_connectivity_new_disk2d()
    else
        throw(ArgumentError("Unknown connectivity"))
    end
    return Connectivity{4}(pointer)
end

function Connectivity{8}(name::Symbol, args...)
    if name == :unitcube
        pointer = p8est_connectivity_new_unitcube()
    elseif name == :periodic
        pointer = p8est_connectivity_new_periodic()
    elseif name == :rotcubes
        pointer = p8est_connectivity_new_rotcubes()
    elseif name == :rotwrap
        pointer = p8est_connectivity_new_rotwrap()
    elseif name == :shell
        pointer = p8est_connectivity_new_shell()
    elseif name == :sphere
        pointer = p8est_connectivity_new_sphere()
    elseif name == :torus
        pointer = p8est_connectivity_new_torus(args...)
    elseif name == :twocubes
        pointer = p8est_connectivity_new_twocubes()
    elseif name == :twotrees
        pointer = p8est_connectivity_new_twotrees(args...)
    elseif name == :twowrap
        pointer = p8est_connectivity_new_twowrap()
    else
        throw(ArgumentError("Unknown connectivity"))
    end
    return Connectivity{8}(pointer)
end


"""
    function Connectivity{X}(vertices::AbstractVector, cells::AbstractVector) where {X}

    Example:
    '''julia
     julia> vertices = [
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 1.0),
                (0.0, 1.0, 0.0),
                (0.0, 1.0, 1.0),
                (1.0, 0.0, 0.0),
                (1.0, 0.0, 1.0),
                (1.0, 1.0, 0.0),
                (1.0, 1.0, 1.0),
            ]
     julia> cells = [
                Int32.((1, 3, 2, 4)),
                Int32.((4, 3, 8, 7)),
                Int32.((7, 5, 8, 6)),
                Int32.((6, 2, 5, 1)),
                Int32.((1, 3, 5, 7)),
                Int32.((2, 4, 6, 8)),
            ]
     julia> conn = Connectivity{4}(vertices, cells)
    '''
       6--------8    cell 1 connects vertices 1, 3, 2, 4
      /        /|    cell 2 connects vertices 4, 3, 8, 7
     /        / |    .
    2--------4  |    .
    |        |  |    .
    |  5     |  7    cell 6 connects vertices 2, 4, 6, 8
    |        | /
    |        |/
    1--------3
"""
function Connectivity{X}(vertices::AbstractVector, cells::AbstractVector) where {X}
    if X == 4
        conn = GC.@preserve vertices cells Connectivity{X}(
            P4est.p4est_connectivity_new(length(vertices), length(cells), 0, 0),
        )
    elseif X == 8
        conn = GC.@preserve vertices cells Connectivity{X}(
            P4est.p8est_connectivity_new(length(vertices), length(cells), 0, 0, 0, 0),
        )
    else
        throw(error("Unsupported cells."))
    end
    trees = unsafe_trees(conn)
    cvertices = unsafe_vertices(conn)
    tree_to_tree = unsafe_tree_to_tree(conn)
    tree_to_face = unsafe_tree_to_face(conn)
    NUM_FACES = (X == 4) ? 4 : 6
    for i in eachindex(cells, trees, tree_to_tree, tree_to_face)
        trees[i] = cells[i] .- 1
        tree_to_tree[i] = ntuple(_ -> (i - 1), NUM_FACES)
        tree_to_face[i] = ntuple(j -> (j - 1), NUM_FACES)
    end

    for i in eachindex(cvertices, vertices)
        cvertices[i] =
            ntuple(j -> (j ≤ length(vertices[i]) ? Float64(vertices[i][j]) : 0.0), Val(3))
    end
    complete!(conn)
    return conn
end

Connectivity{4}(name::String) = Connectivity{4}(p4est_connectivity_read_inp(name))
Connectivity{8}(name::String) = Connectivity{8}(p8est_connectivity_read_inp(name))

function Base.unsafe_convert(
    ::Type{Ptr{p4est_connectivity}},
    c::Connectivity{4,Ptr{p4est_connectivity}},
)
    return c.pointer
end
function Base.unsafe_convert(
    ::Type{Ptr{p8est_connectivity}},
    c::Connectivity{8,Ptr{p8est_connectivity}},
)
    return c.pointer
end

function Base.:(==)(x::Connectivity{4}, y::Connectivity{4})
    p4est_connectivity_is_equivalent(x.pointer, y.pointer) == 1
end
function Base.:(==)(x::Connectivity{8}, y::Connectivity{8})
    p8est_connectivity_is_equivalent(x.pointer, y.pointer) == 1
end

Base.isvalid(c::Connectivity{4}) =
    GC.@preserve c p4est_connectivity_is_valid(c.pointer) == 1
Base.sizeof(c::Connectivity{4}) = Int(p4est_connectivity_memory_used(c.pointer))

Base.isvalid(c::Connectivity{8}) = p8est_connectivity_is_valid(c.pointer) == 1
Base.sizeof(c::Connectivity{8}) = Int(p8est_connectivity_memory_used(c.pointer))

function unsafe_vertices(c::Connectivity{X}) where {X}
    cs = PointerWrapper(c.pointer)
    v = unsafe_wrap(
        Matrix{Cdouble},
        pointer(cs.vertices),
        (3, Int(cs.num_vertices[])),
        own = false,
    )
    return reinterpret(reshape, NTuple{3,Cdouble}, v)
end

function unsafe_trees(c::Connectivity{X}) where {X}
    cs = PointerWrapper(c.pointer)
    ttv = unsafe_wrap(
        Matrix{p4est_topidx_t},
        pointer(cs.tree_to_vertex),
        (X, Int(cs.num_trees[])),
        own = false,
    )
    return reinterpret(reshape, NTuple{X,p4est_topidx_t}, ttv)
end

function unsafe_tree_to_tree(c::Connectivity{X}) where {X}
    cs = PointerWrapper(c.pointer)
    sides = X == 4 ? 4 : 6
    ttt = unsafe_wrap(
        Matrix{p4est_topidx_t},
        pointer(cs.tree_to_tree),
        (sides, Int(cs.num_trees[])),
        own = false,
    )
    return reinterpret(reshape, NTuple{sides,p4est_topidx_t}, ttt)
end

function unsafe_tree_to_face(c::Connectivity{X}) where {X}
    cs = PointerWrapper(c.pointer)
    sides = X == 4 ? 4 : 6
    ttf = unsafe_wrap(
        Matrix{Int8},
        pointer(cs.tree_to_face),
        (sides, Int(cs.num_trees[])),
        own = false,
    )
    return reinterpret(reshape, NTuple{sides,Int8}, ttf)
end

function unsafe_tree_to_corner(c::Connectivity{X}) where {X}
    cs = PointerWrapper(c.pointer)
    s = pointer(cs.tree_to_corner) == C_NULL ? 0 : Int(cs.num_trees[])
    ttc =
        unsafe_wrap(Matrix{p4est_topidx_t}, pointer(cs.tree_to_corner), (X, s), own = false)
    return reinterpret(reshape, NTuple{X,p4est_topidx_t}, ttc)
end

function unsafe_ctt_offset(c::Connectivity{X}) where {X}
    cs = PointerWrapper(c.pointer)
    ctt = unsafe_wrap(
        Vector{p4est_topidx_t},
        pointer(cs.ctt_offset),
        (Int(cs.num_corners[] + 1),),
        own = false,
    )
    return ctt
end

function unsafe_corner_to_tree_array(c::Connectivity{X}) where {X}
    cs = PointerWrapper(c.pointer)
    off = unsafe_ctt_offset(c)
    s = length(off) > 0 ? last(off) : 0
    return unsafe_wrap(
        Vector{p4est_topidx_t},
        pointer(cs.corner_to_tree),
        (s,),
        own = false,
    )
end

function unsafe_corner_to_corner_array(c::Connectivity{X}) where {X}
    cs = PointerWrapper(c.pointer)
    off = unsafe_ctt_offset(c)
    s = length(off) > 0 ? last(off) : 0
    return unsafe_wrap(Vector{Int8}, pointer(cs.corner_to_corner), (s,), own = false)
end

function unsafe_corner_to_tree(c::Connectivity{X}) where {X}
    cs = PointerWrapper(c.pointer)
    o = unsafe_ctt_offset(c)
    ctt = P4estTypes.unsafe_corner_to_tree_array(c)
    ctc = P4estTypes.unsafe_corner_to_corner_array(c)
    s = SparseMatrixCSC(cs.num_trees[], cs.num_corners[], o .+ 1, ctt .+ 1, ctc)
    return s
end

"""
    refine(c::Connectivity{4}, nedge)

Returns a new [`Connectivity`](@ref) that is `c` uniformly refined with `nedge`
new trees in each direction.
"""
function refine(c::Connectivity{4}, nedge)
    return GC.@preserve c nedge Connectivity{4}(p4est_connectivity_refine(c.pointer, nedge))
end

"""
    refine(c::Connectivity{8}, nedge)

Returns a new [`Connectivity`](@ref) that is `c` uniformly refined with `nedge`
new trees in each direction.
"""
function refine(c::Connectivity{8}, nedge)
    return Connectivity{8}(p8est_connectivity_refine(c.pointer, nedge))
end

reduce!(c::Connectivity{4}) = GC.@preserve c p4est_connectivity_reduce(c.pointer)
reduce!(c::Connectivity{8}) = p8est_connectivity_reduce(c.pointer)

complete!(c::Connectivity{4}) = GC.@preserve c p4est_connectivity_complete(c.pointer)
complete!(c::Connectivity{8}) = p8est_connectivity_complete(c.pointer)

"""
    brick(n::NTuple{2, Integer}, p::NTuple{2, Bool}=(false, false))

Returns a new [`Connectivity`](@ref) that is a rectangular `n[1]`-by-`n[2]`
quadtree connectivity.  The brick is periodic in x and y if `p[1]` and `p[2]`
are `true`, respectively.
"""
function brick(n::Tuple{Integer,Integer}, p::Tuple{Bool,Bool} = (false, false))
    return Connectivity{4}(p4est_connectivity_new_brick(n..., p...))
end

"""
    brick(n::NTuple{3, Integer}, p::NTuple{3, Bool}=(false, false, false))

Returns a new [`Connectivity`](@ref) that is a rectangular
`n[1]`-by-`n[2]`-by-`n[3]` octree mesh.  The brick is periodic in x, y, and z
if `p[1]`, `p[2]`, and `p[3]` are `true`, respectively.
"""
function brick(
    n::Tuple{Integer,Integer,Integer},
    p::Tuple{Bool,Bool,Bool} = (false, false, false),
)
    return Connectivity{8}(p8est_connectivity_new_brick(n..., p...))
end

"""
    brick(l, m, p=false, q=false)

Returns a new [`Connectivity`](@ref) that is a rectangular `l`-by-`m` quadtree
mesh.  The brick is periodic in x and y if `p` and `q` are `true`, respectively.
"""
brick(l::Integer, m::Integer, p::Bool = false, q::Bool = false) = brick((l, m), (p, q))

"""
    brick(l, m, n, p=false, q=false, r=false)

Returns a new [`Connectivity`](@ref) that is a rectangular `l`-by-`m`-by-`n`
octree mesh.  The brick is periodic in x, y, and z if `p`, `q`, and `r` are
`true`, respectively.
"""
function brick(
    l::Integer,
    m::Integer,
    n::Integer,
    p::Bool = false,
    q::Bool = false,
    r::Bool = false,
)
    return brick((l, m, n), (p, q, r))
end

function Base.show(io::IO, ::Connectivity{X}) where {X}
    print(io, "Connectivity{", string(X), "}")
end

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, c::Connectivity{X}) where {X}
    show(io, c)
    if !get(io, :compact, false)
        GC.@preserve c begin
            # Should we print everything 1 based?
            println(io, "\nnote: the following entries are zero-based")
            println(io, "\ntrees:")
            show(io, mime, unsafe_trees(c))
            println(io, "\nvertices:")
            show(io, mime, unsafe_vertices(c))
            println(io, "\ntree to tree:")
            show(io, mime, unsafe_tree_to_tree(c))
            println(io, "\ntree to face:")
            show(io, mime, unsafe_tree_to_face(c))
            println(io, "\ntree to corner:")
            show(io, mime, unsafe_tree_to_corner(c))
            println(io, "\ncorners:")
            show(io, mime, unsafe_corner_to_tree(c))
        end
    end
end

@recipe function f(c::Connectivity{4})
    verts = unsafe_vertices(c)

    xlims = extrema(getindex.(verts, 1))
    ylims = extrema(getindex.(verts, 2))
    zlims = extrema(getindex.(verts, 3))

    isconstz = zlims[1] == zlims[2]

    xlabel --> "x"
    ylabel --> "y"
    zlabel --> "z"

    aspect_ratio --> :equal
    legend --> false
    grid --> false

    @series begin
        seriestype --> :path
        linecolor --> :gray
        linewidth --> 1

        x = []
        y = []
        z = []
        for tree in unsafe_trees(c)
            for i in (1, 2, 4, 3, 1)
                xi, yi, zi = verts[tree[i]+1]

                push!(x, xi)
                push!(y, yi)
                push!(z, zi)
            end

            push!(x, NaN)
            push!(y, NaN)
            push!(z, NaN)
        end
        if isconstz
            x, y
        else
            x, y, z
        end
    end
end

@recipe function f(c::Connectivity{8})
    verts = unsafe_vertices(c)

    xlims = extrema(getindex.(verts, 1))
    ylims = extrema(getindex.(verts, 2))
    zlims = extrema(getindex.(verts, 3))

    xlabel --> "x"
    ylabel --> "y"
    zlabel --> "z"

    aspect_ratio --> :equal
    legend --> false
    grid --> false

    @series begin
        seriestype --> :path
        linecolor --> :gray
        linewidth --> 1

        x = []
        y = []
        z = []
        for tree in unsafe_trees(c)
            for j in (0, 4)
                for i in (1 + j, 2 + j, 4 + j, 3 + j, 1 + j)
                    xi, yi, zi = verts[tree[i]+1]

                    push!(x, xi)
                    push!(y, yi)
                    push!(z, zi)
                end

                push!(x, NaN)
                push!(y, NaN)
                push!(z, NaN)
            end

            for j = 0:3
                for i in (1 + j, 5 + j)
                    xi, yi, zi = verts[tree[i]+1]

                    push!(x, xi)
                    push!(y, yi)
                    push!(z, zi)
                end

                push!(x, NaN)
                push!(y, NaN)
                push!(z, NaN)
            end
        end

        x, y, z
    end
end
