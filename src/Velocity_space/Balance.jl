"""
$(TYPEDSIGNATURES)
Enforce the 2:1 balance constraint in velocity space for two neighboring physical cells.

Walks both velocity-space grids simultaneously (in their shared DFS tree order) and
refines any local cell whose level is more than 1 below the corresponding neighbor cell's
level.  The neighbor may be a local [`VsData`](@ref) or a read-only
[`GhostVsData`](@ref); only the local grid is ever modified.

Returns `true` if at least one local cell was refined.
"""
function vs_balance_pair!(
    vs_local::VsData{DIM,NDF},
    vs_neighbor::AbstractVsData{DIM,NDF},
    ds::AbstractVector,
) where {DIM,NDF}
    changed = false
    refine_flags = Bool[]
    # Refine local cell i one level whenever the maximum neighbor level over its spatial extent
    # exceeds level(i)+1 (the 2:1 violation).  A read-only two-pointer walk fills the per-cell
    # flags; a single streaming refine then applies them.  Looping to a fixed point fully
    # balances this pair (a cell may need several levels), each level costing O(N) instead of the
    # O(N) positional insert! per refined cell.
    while true
        nL = vs_local.vs_num;
        nN = vs_neighbor.vs_num
        resize!(refine_flags, nL);
        fill!(refine_flags, false)
        i = 1;
        j = 1
        flag = 0.0          # fraction of the current (coarser) neighbor cell already covered
        @inbounds while i <= nL && j <= nN
            l = vs_local.level[i];
            ln = vs_neighbor.level[j]
            if l == ln
                i += 1;
                j += 1
            elseif l < ln
                # Local cell i is coarser: accumulate the neighbor cells covering it, tracking
                # their maximum level, then decide.
                maxlev = 0;
                f = 0.0
                while f < 1.0 - 1e-10 && j <= nN
                    lv = vs_neighbor.level[j]
                    lv > maxlev && (maxlev = Int(lv))
                    f += 1.0 / 2^(DIM * (lv - l))
                    j += 1
                end
                refine_flags[i] = maxlev > l + 1
                i += 1
            else
                # Local cell i is finer: one neighbor cell covers it (and more); no violation.
                flag += 1.0 / 2^(DIM * (l - ln))
                if flag >= 1.0 - 1e-10
                    j += 1;
                    flag = 0.0
                end
                i += 1
            end
        end
        any(refine_flags) || break
        refine_grid_stream!(vs_local, refine_flags, ds)
        changed = true
    end
    return changed
end

"""
$(TYPEDSIGNATURES)
Single-rank velocity-space balance pass.

Iterates over all local physical cells and their neighbors (local *and* ghost) until no
further violations remain among local cells.  Returns `true` if at least one local
velocity cell was refined.

Ghost cells are read-only: violations where a ghost is coarser than the local cell are
resolved on the remote rank in subsequent communication rounds.
"""
function vs_balance_local!(ka::KA{DIM,NDF}) where {DIM,NDF}
    !isa(ka.kinfo.config.quadrature, Vector) && return false

    trees = ka.kdata.field.trees
    ds = [
        (ka.kinfo.config.quadrature[2*i] - ka.kinfo.config.quadrature[2*i-1]) /
        ka.kinfo.config.vs_trees_num[i] for i = 1:DIM
    ]

    any_changed = false
    changed = true
    while changed
        changed = false
        for tree in trees.data
            for ps_data in tree
                isa(ps_data, InsideSolidData) && continue
                ps_data.bound_enc < 0 && continue

                for face_idx = 1:(2*DIM)
                    for neighbor in ps_data.neighbor.data[face_idx]
                        neighbor === nothing && continue
                        isa(neighbor, AbstractInsideSolidData) && continue
                        isa(neighbor, SolidNeighbor) && continue

                        if vs_balance_pair!(ps_data.vs_data, neighbor.vs_data, ds)
                            changed = true
                            any_changed = true
                        end
                    end
                end
            end
        end
    end

    return any_changed
end

"""
$(TYPEDSIGNATURES)
After `ghost_wrap` has been rebuilt, patch the stale ghost-cell pointers stored in
every local cell's `neighbor.data`.

`old_ghost_wrap` and `new_ghost_wrap` must have the same length and the same ordering
(guaranteed when the p4est ghost topology has not changed).  The patch is done by
position: `old_ghost_wrap[i]` is replaced by `new_ghost_wrap[i]`.
"""
function _patch_neighbor_ghost_refs!(
    ka::KA{DIM,NDF},
    old_ghost_wrap::AbstractVector,
    new_ghost_wrap::AbstractVector,
) where {DIM,NDF}
    # Map Julia object identity (old ghost) → new ghost object.
    ghost_patch = Dict{UInt64,AbstractGhostPsData{DIM,NDF}}()
    sizehint!(ghost_patch, length(old_ghost_wrap))
    for i in eachindex(old_ghost_wrap)
        ghost_patch[objectid(old_ghost_wrap[i])] = new_ghost_wrap[i]
    end

    for tree in ka.kdata.field.trees.data
        for ps_data in tree
            isa(ps_data, InsideSolidData) && continue
            for face_idx in eachindex(ps_data.neighbor.data)
                face_neighbors = ps_data.neighbor.data[face_idx]
                for k in eachindex(face_neighbors)
                    nb = face_neighbors[k]
                    nb === nothing && continue
                    isa(nb, AbstractGhostPsData) || continue
                    new_nb = get(ghost_patch, objectid(nb), nothing)
                    new_nb === nothing && continue
                    face_neighbors[k] = new_nb
                end
            end
        end
    end
    return nothing
end

"""
$(TYPEDSIGNATURES)
Re-communicate the velocity-space structure (levels, midpoints, weights, and
distribution-function values) of ghost cells after local velocity-space refinement,
without rebuilding the physical-mesh ghost topology.

Compared to `update_ghost!`, this function skips the `p4est_ghost_destroy / p4est_ghost_new` calls because the physical quadrant adjacency is unchanged.  Ghost-cell
pointers cached in every local cell's `neighbor.data` are patched in-place so that they
refer to the freshly allocated ghost objects.

The p4est mesh and the face list are **not** rebuilt here; call `update_neighbor!` and
`update_faces!` (or the combined `amr_recover!`) after the full balance loop if those
structures must be up-to-date for subsequent operations.
"""
function vs_ghost_exchange!(p4est::P_pxest_t, ka::KA{DIM,NDF}) where {DIM,NDF}
    kinfo = ka.kinfo
    # Hold a reference to the old ghost_wrap so that its objects stay alive (and their
    # Julia objectid remains valid) until after _patch_neighbor_ghost_refs! completes.
    old_ghost_wrap = ka.kdata.ghost.ghost_wrap
    # Re-exchange velocity-space structure; Julia GC handles old buffer cleanup.
    new_gb, new_gi = initialize_ghost_pool(p4est, kinfo)
    ka.kdata.ghost.ghost_buffer = new_gb
    ka.kdata.ghost.ghost_info = new_gi
    ka.kdata.ghost.ghost_wrap = initialize_ghost_wrap(kinfo, new_gb, new_gi)
    # Replace stale ghost references in neighbor.data with the new objects.
    _patch_neighbor_ghost_refs!(ka, old_ghost_wrap, ka.kdata.ghost.ghost_wrap)
    return nothing
end

"""
$(TYPEDSIGNATURES)
Enforce the global 2:1 balance constraint in velocity space across all neighboring
physical cells, including pairs that span different MPI ranks.

**Algorithm**

Each iteration:

 1. Runs a local balance pass ([`vs_balance!`](@ref)) that refines local velocity cells
    which are more than 1 level coarser than any neighbor (local or ghost).
 2. Performs an `MPI_Allreduce` to check whether *any* rank made changes.
 3. If changes occurred on at least one rank, calls [`vs_ghost_exchange!`](@ref) so that
    every rank learns the updated velocity-space levels of its neighbors before the next
    pass.

The loop terminates when no rank reports a change, guaranteeing that the 2:1 constraint
is satisfied globally.

For single-rank runs the MPI steps are skipped.

**Note** After this function returns:

  - `ghost_wrap` and `neighbor.data` contain up-to-date velocity-space level data.
  - The p4est mesh and the face list are **not** rebuilt.  Call `update_neighbor!(p4est, ka)`
    and `update_faces!(p4est, ka)` (or `amr_recover!(p4est, ka)`) afterwards if those
    structures are needed for flux or slope calculations.
"""
function vs_balance!(ka::KA{DIM,NDF}) where {DIM,NDF}
    p4est = ka.kinfo.forest.p4est
    !isa(ka.kinfo.config.quadrature, Vector) && return false

    if MPI.Comm_size(MPI.COMM_WORLD) == 1
        return vs_balance_local!(ka)
    end

    any_changed = false
    changed_global = true
    while changed_global
        changed_local = vs_balance_local!(ka)
        # Reduce across all ranks: any change on any rank keeps the loop going.
        changed_global = Bool(MPI.Allreduce(Int(changed_local), +, MPI.COMM_WORLD) > 0)
        any_changed |= changed_global
        # Rebuild ghost velocity-space data so that each rank sees its neighbours'
        # updated levels in the next pass.
        changed_global && vs_ghost_exchange!(p4est, ka)
    end

    return any_changed
end
