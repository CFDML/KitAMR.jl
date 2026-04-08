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
    lnmidpoint = reshape(vs_local.midpoint, :)
    lndf       = reshape(vs_local.df, :)
    lnsdf      = reshape(vs_local.sdf, :)
    lnflux     = reshape(vs_local.flux, :)

    changed = false
    i = 1          # current index in local grid
    j = 1          # current index in neighbor grid
    flag = 0.0     # fraction of the current neighbor cell already accounted for

    midpoint_index = Vector{Int}(undef, DIM)
    df_index       = Vector{Int}(undef, NDF)

    while i <= vs_local.vs_num && j <= vs_neighbor.vs_num
        l  = vs_local.level[i]
        ln = vs_neighbor.level[j]

        if l == ln
            # Same level – one-to-one correspondence, no violation.
            i += 1
            j += 1

        elseif l < ln
            # Local cell i is coarser than the corresponding neighbor region.
            if ln - l > 1
                # Violation: refine the local cell one level.
                for k in 1:DIM
                    midpoint_index[k] = (k - 1) * vs_local.vs_num + i
                end
                for k in 1:NDF
                    df_index[k] = (k - 1) * vs_local.vs_num + i
                end
                midpoint_i   = lnmidpoint[midpoint_index]
                df_i         = lndf[df_index]
                midpoint_new = midpoint_refine(DIM, midpoint_i, l, ds)
                df_new       = df_refine(DIM, midpoint_i, midpoint_new, df_i)
                vs_local.vs_num += 2^DIM - 1
                level_refine_replace!(DIM, vs_local.level, i)
                weight_refine_replace!(DIM, vs_local.weight, i)
                midpoint_refine_replace!(DIM, lnmidpoint, midpoint_new, vs_local.vs_num, i)
                df_refine_replace!(DIM, NDF, lndf, df_new, vs_local.vs_num, i)
                sdf_refine_replace!(DIM, NDF, lnsdf)
                flux_refine_replace!(DIM, NDF, lnflux)
                changed = true
                # Do NOT advance i: re-examine the newly created children against j.
            else
                # ln == l + 1: currently acceptable, but neighbor sub-cells may be
                # refined further.  Check every neighbor sub-cell before consuming it;
                # if any has level > l + 1, refine local cell i first and re-examine.
                local_violated = false
                j_scan = j
                flag_scan = flag
                while flag_scan < 1.0 - 1e-10 && j_scan <= vs_neighbor.vs_num
                    if vs_neighbor.level[j_scan] - l > 1
                        local_violated = true
                        break
                    end
                    flag_scan += 1.0 / 2^(DIM * (vs_neighbor.level[j_scan] - l))
                    j_scan += 1
                end
                if local_violated
                    # Refine local cell i and re-examine without advancing j.
                    for k in 1:DIM
                        midpoint_index[k] = (k - 1) * vs_local.vs_num + i
                    end
                    for k in 1:NDF
                        df_index[k] = (k - 1) * vs_local.vs_num + i
                    end
                    midpoint_i   = lnmidpoint[midpoint_index]
                    df_i         = lndf[df_index]
                    midpoint_new = midpoint_refine(DIM, midpoint_i, l, ds)
                    df_new       = df_refine(DIM, midpoint_i, midpoint_new, df_i)
                    vs_local.vs_num += 2^DIM - 1
                    level_refine_replace!(DIM, vs_local.level, i)
                    weight_refine_replace!(DIM, vs_local.weight, i)
                    midpoint_refine_replace!(DIM, lnmidpoint, midpoint_new, vs_local.vs_num, i)
                    df_refine_replace!(DIM, NDF, lndf, df_new, vs_local.vs_num, i)
                    sdf_refine_replace!(DIM, NDF, lnsdf)
                    flux_refine_replace!(DIM, NDF, lnflux)
                    changed = true
                    # Do NOT advance i or j: re-examine the new children against j.
                else
                    # All neighbor sub-cells satisfy the constraint; consume them.
                    while flag < 1.0 - 1e-10 && j <= vs_neighbor.vs_num
                        flag += 1.0 / 2^(DIM * (vs_neighbor.level[j] - l))
                        j += 1
                    end
                    flag = 0.0
                    i += 1
                end
            end

        else
            # l > ln: local cell is finer than the neighbor region – no local refinement
            # needed.  Track what fraction of neighbor cell j has been covered.
            flag += 1.0 / 2^(DIM * (l - ln))
            if flag >= 1.0 - 1e-10
                j += 1
                flag = 0.0
            end
            i += 1
        end
    end

    vs_local.midpoint = reshape(lnmidpoint, vs_local.vs_num, DIM)
    vs_local.df       = reshape(lndf,       vs_local.vs_num, NDF)
    vs_local.sdf      = reshape(lnsdf,      vs_local.vs_num, NDF, DIM)
    vs_local.flux     = reshape(lnflux,     vs_local.vs_num, NDF)

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
        ka.kinfo.config.vs_trees_num[i]
        for i in 1:DIM
    ]

    any_changed = false
    changed = true
    while changed
        changed = false
        for tree in trees.data
            for ps_data in tree
                isa(ps_data, InsideSolidData) && continue
                ps_data.bound_enc < 0 && continue

                for face_idx in 1:2*DIM
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
    ghost_patch = Dict{UInt64, AbstractGhostPsData{DIM,NDF}}()
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

Compared to `update_ghost!`, this function skips the `p4est_ghost_destroy /
p4est_ghost_new` calls because the physical quadrant adjacency is unchanged.  Ghost-cell
pointers cached in every local cell's `neighbor.data` are patched in-place so that they
refer to the freshly allocated ghost objects.

The p4est mesh and the face list are **not** rebuilt here; call `update_neighbor!` and
`update_faces!` (or the combined `amr_recover!`) after the full balance loop if those
structures must be up-to-date for subsequent operations.
"""
function vs_ghost_exchange!(p4est::P_pxest_t, ka::KA{DIM,NDF}) where {DIM,NDF}
    kinfo = ka.kinfo
    # Hold a reference to the old ghost_wrap so that its objects stay alive (and their
    # Julia objectid remains valid) after the C-level buffers are freed.
    old_ghost_wrap = ka.kdata.ghost.ghost_wrap
    # Free old C-level buffers (ghost_datas, ghost_slopes, ghost_structures, mirrors).
    finalize_ghost!(ka.kdata.ghost.ghost_pointers)
    # Re-exchange: get_mirror_data internally calls get_vs_num which updates
    # kinfo.status.max_vs_num via MPI_Allreduce, handling any increase in vs_num.
    ka.kdata.ghost.ghost_pointers = initialize_ghost_pointers(p4est, kinfo)
    ka.kdata.ghost.ghost_wrap     = initialize_ghost_wrap(kinfo, ka.kdata.ghost.ghost_pointers)
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
    !isa(ka.kinfo.config.quadrature, Vector) && return nothing

    if MPI.Comm_size(MPI.COMM_WORLD) == 1
        vs_balance_local!(ka)
        return nothing
    end

    changed_global = true
    while changed_global
        if MPI.Comm_rank(MPI.COMM_WORLD) == 0
            @show "vs balance communication..."
        end
        changed_local  = vs_balance_local!(ka)
        # Reduce across all ranks: any change on any rank keeps the loop going.
        changed_global = Bool(MPI.Allreduce(Int(changed_local), +, MPI.COMM_WORLD) > 0)
        # Rebuild ghost velocity-space data so that each rank sees its neighbours'
        # updated levels in the next pass.
        changed_global && vs_ghost_exchange!(p4est, ka)
    end

    return nothing
end
