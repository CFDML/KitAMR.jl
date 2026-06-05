# Minimal local reimplementation of the small subset of `P4estTypes.jl` that KitAMR
# actually uses (only in `Mesh/Connectivity.jl`): the `Connectivity` wrapper, the
# `(vertices, cells)` constructor, `brick`, and `unsafe_vertices`. This removes the
# dependency on the vendored `lib/P4estTypes`. The underlying C structs/functions come
# from the `P4est` bindings, which are already brought in via `using P4est`.

"""
$(TYPEDEF)
Lightweight wrapper around a p4est (`X == 4`, 2D) or p8est (`X == 8`, 3D) connectivity
C struct. `pointer` is a `Ptr{p4est_connectivity}` / `Ptr{p8est_connectivity}` suitable
for passing to `p4est_new_ext` / `p8est_new_ext`.

No finalizer is attached: the connectivity must outlive the forest built from it and is
released when the process exits, matching the upstream `P4estTypes` behavior.
"""
mutable struct Connectivity{X,P}
    pointer::P
    Connectivity{4}(pointer::Ptr{p4est_connectivity}) = new{4,typeof(pointer)}(pointer)
    Connectivity{8}(pointer::Ptr{p8est_connectivity}) = new{8,typeof(pointer)}(pointer)
end

# --- writable, non-owning views into the connectivity arrays -----------------------
"""
$(TYPEDSIGNATURES)
Return a writable view (length `num_vertices`) of the connectivity's vertex coordinates
as `NTuple{3,Cdouble}`. Writing through it mutates the underlying C memory.
"""
function unsafe_vertices(c::Connectivity)
    cs = PointerWrapper(c.pointer)
    v = unsafe_wrap(Matrix{Cdouble}, pointer(cs.vertices), (3, Int(cs.num_vertices[])), own = false)
    return reinterpret(reshape, NTuple{3,Cdouble}, v)
end
function unsafe_trees(c::Connectivity{X}) where {X}
    cs = PointerWrapper(c.pointer)
    ttv = unsafe_wrap(Matrix{p4est_topidx_t}, pointer(cs.tree_to_vertex), (X, Int(cs.num_trees[])), own = false)
    return reinterpret(reshape, NTuple{X,p4est_topidx_t}, ttv)
end
function unsafe_tree_to_tree(c::Connectivity{X}) where {X}
    cs = PointerWrapper(c.pointer)
    sides = X == 4 ? 4 : 6
    ttt = unsafe_wrap(Matrix{p4est_topidx_t}, pointer(cs.tree_to_tree), (sides, Int(cs.num_trees[])), own = false)
    return reinterpret(reshape, NTuple{sides,p4est_topidx_t}, ttt)
end
function unsafe_tree_to_face(c::Connectivity{X}) where {X}
    cs = PointerWrapper(c.pointer)
    sides = X == 4 ? 4 : 6
    ttf = unsafe_wrap(Matrix{Int8}, pointer(cs.tree_to_face), (sides, Int(cs.num_trees[])), own = false)
    return reinterpret(reshape, NTuple{sides,Int8}, ttf)
end

complete!(c::Connectivity{4}) = GC.@preserve c p4est_connectivity_complete(c.pointer)
complete!(c::Connectivity{8}) = GC.@preserve c p8est_connectivity_complete(c.pointer)

# --- build a connectivity from a vertex list and 1-based cell -> vertex tuples ------
"""
$(TYPEDSIGNATURES)
Build a [`Connectivity`](@ref) from `vertices` (each an `NTuple` of coordinates) and
`cells` (each a 1-based `NTuple{X}` of vertex indices). `X == 4` for 2D quads, `X == 8`
for 3D hexes. Self-connecting faces are filled in and `complete!` resolves the topology.
"""
function Connectivity{X}(vertices::AbstractVector, cells::AbstractVector) where {X}
    if X == 4
        conn = GC.@preserve vertices cells Connectivity{X}(
            p4est_connectivity_new(length(vertices), length(cells), 0, 0))
    elseif X == 8
        conn = GC.@preserve vertices cells Connectivity{X}(
            p8est_connectivity_new(length(vertices), length(cells), 0, 0, 0, 0))
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
        cvertices[i] = ntuple(j -> (j ≤ length(vertices[i]) ? Float64(vertices[i][j]) : 0.0), Val(3))
    end
    complete!(conn)
    return conn
end

# --- rectangular brick connectivities ----------------------------------------------
"""
$(TYPEDSIGNATURES)
Rectangular `n[1]`-by-`n[2]` quadtree (2D) connectivity, periodic in x/y where `p` is true.
"""
brick(n::Tuple{Integer,Integer}, p::Tuple{Bool,Bool} = (false, false)) =
    Connectivity{4}(p4est_connectivity_new_brick(n..., p...))
"""
$(TYPEDSIGNATURES)
Rectangular `n[1]`-by-`n[2]`-by-`n[3]` octree (3D) connectivity, periodic in x/y/z where `p` is true.
"""
brick(n::Tuple{Integer,Integer,Integer}, p::Tuple{Bool,Bool,Bool} = (false, false, false)) =
    Connectivity{8}(p8est_connectivity_new_brick(n..., p...))



function Cartesian_connectivity(Nx, Ny, xmin, xmax, ymin, ymax)
    vertices_C = Array{NTuple{2,Float64}}(undef, Nx + 1, Ny + 1)
    dx = (xmax - xmin) / Nx
    dy = (ymax - ymin) / Ny
    for j = 1:Ny+1
        for i = 1:Nx+1
            vertices_C[i, j] = (xmin + (i - 1) * dx, ymin + (j - 1) * dy)
        end
    end
    vL = LinearIndices(vertices_C)
    vertices = reshape(vertices_C, :)
    cell2ver = Array{NTuple{4,Int32}}(undef, Nx, Ny)
    for j in axes(cell2ver, 2)
        for i in axes(cell2ver, 1)
            cell2ver[i, j] = (vL[i, j], vL[i+1, j], vL[i, j+1], vL[i+1, j+1])
        end
    end
    cells = reshape(cell2ver, :)
    connectivity = GC.@preserve vertices cells Connectivity{4}(vertices, cells)
end
function Cartesian_connectivity(Nx,Ny,Nz,xmin,xmax,ymin,ymax,zmin,zmax)
    vertices_C = Array{NTuple{3,Float64}}(undef, Nx + 1, Ny + 1, Nz + 1)
    dx = (xmax - xmin) / Nx
    dy = (ymax - ymin) / Ny
    dz = (zmax - zmin) / Nz
    for k = 1:Nz+1
        for j = 1:Ny+1
            for i = 1:Nx+1
                vertices_C[i, j, k] = (xmin + (i - 1) * dx, ymin + (j - 1) * dy, zmin + (k - 1) * dz)
            end
        end
    end
    vL = LinearIndices(vertices_C)
    vertices = reshape(vertices_C, :)
    cell2ver = Array{NTuple{8,Int32}}(undef, Nx, Ny, Nz)
    for k in axes(cell2ver, 3)
        for j in axes(cell2ver, 2)
            for i in axes(cell2ver, 1)
                cell2ver[i, j, k] = (vL[i, j, k], vL[i+1, j, k], vL[i, j+1, k], vL[i+1, j+1, k],
                                     vL[i, j, k+1], vL[i+1, j, k+1], vL[i, j+1, k+1], vL[i+1, j+1, k+1])
            end
        end
    end
    cells = reshape(cell2ver, :)
    connectivity = GC.@preserve vertices cells Connectivity{8}(vertices, cells)
end
function set_connectivity(kinfo::KInfo{DIM}) where{DIM}
    domain = kinfo.config.domain

    periodic_dirs = ntuple(i -> begin
        first_boundary = (i-1)*2 + 1  # 1, 3, 5 for i = 1, 2, 3
        second_boundary = first_boundary + 1  # 3, 5, 7 for i = 1, 2, 3
        nameof(typeof(domain[first_boundary]).parameters[1]) == :Period &&
        nameof(typeof(domain[second_boundary]).parameters[1]) == :Period
    end, DIM)
    
    if periodic_dirs != ntuple(_ -> false, DIM)
        connectivity_ps = set_periodic_connectivity(kinfo, periodic_dirs)
    else
        connectivity_ps = Cartesian_connectivity(kinfo.config.trees_num..., kinfo.config.geometry...)
    end
    return connectivity_ps
end

function set_periodic_connectivity(kinfo::KInfo{DIM}, periodic_dirs) where{DIM}
    geometry = kinfo.config.geometry
    bounds = if DIM == 2
        ntuple(i -> begin
            if i <= 2
                min_idx = (i-1)*2 + 1  # 1, 3 for i = 1, 2
                max_idx = min_idx + 1   # 2, 4 for i = 1, 2
                (geometry[min_idx], geometry[max_idx])  # (xmin,xmax), (ymin,ymax)
            else
                (0.0, 0.0)
            end
        end, 3)
    else
        ntuple(i -> begin
            min_idx = (i-1)*2 + 1  # 1, 3, 5 for i = 1, 2, 3
            max_idx = min_idx + 1   # 2, 4, 6 for i = 1, 2, 3
            (geometry[min_idx], geometry[max_idx])  # (xmin,xmax), (ymin,ymax), (zmin,zmax)
        end, 3)
    end
    
    connectivity_ps = brick(
        Tuple(kinfo.config.trees_num),
        periodic_dirs
    )
    vertices = unsafe_vertices(connectivity_ps)
    for i in eachindex(vertices)
        normalized_coords = vertices[i]
        
        vertices[i] = ntuple(j -> begin
            if j <= DIM
                normalized_coords[j] * (bounds[j][2] - bounds[j][1])/kinfo.config.trees_num[j] + bounds[j][1]
            else
                0.0
            end
        end, 3)
    end
    return connectivity_ps
end