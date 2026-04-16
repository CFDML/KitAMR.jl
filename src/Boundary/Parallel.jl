"""
Sparse variable-size MPI exchange where `ghost_elem_sizes[i] == 0` marks cells to skip.
Unlike `_mpi_exchange_varsize!`, the active ghost cells from a given rank are NOT
contiguous in `ghost_buffer`, so receives are packed into a temporary buffer per source
rank and then scattered into the correct positions after `MPI.Waitall`.
"""
function _mpi_solid_exchange_varsize!(
    ghost::P_pxest_ghost_t,
    ghost_buffer::Vector{Float64},
    ghost_offsets::Vector{Int},
    ghost_elem_sizes::Vector{Int},   # 0 for non-solid ghosts
    mirror_data_bufs::Vector{Vector{Float64}},
    mirror_elem_sizes::Vector{Int},  # 0 for non-solid mirrors
    tag::Int,
)
    mpisize, proc_offsets, mirror_proc_offsets, mirror_proc_mirrors =
        _ghost_comm_arrays(ghost)
    myrank = MPI.Comm_rank(MPI.COMM_WORLD)
    comm   = MPI.COMM_WORLD

    reqs      = Vector{MPI.Request}(undef, 0)
    send_bufs = Vector{Vector{Float64}}(undef, 0)
    # Each entry: (recv_buf, solid_ghost_indices, solid_ghost_sizes) for one source rank.
    recv_infos = Vector{Tuple{Vector{Float64}, Vector{Int}, Vector{Int}}}(undef, 0)

    # --- Post receives grouped by source rank (only solid ghost cells) ---
    for r in 0:mpisize-1
        r == myrank && continue
        g_start = Int(proc_offsets[r+1]) + 1
        g_end   = Int(proc_offsets[r+2])
        g_start > g_end && continue
        solid_idx  = Int[]
        solid_szs  = Int[]
        recv_elems = 0
        for i in g_start:g_end
            sz = ghost_elem_sizes[i]
            sz == 0 && continue
            push!(solid_idx, i)
            push!(solid_szs, sz)
            recv_elems += sz
        end
        recv_elems == 0 && continue
        recv_buf = Vector{Float64}(undef, recv_elems)
        push!(recv_infos, (recv_buf, solid_idx, solid_szs))
        push!(reqs, MPI.Irecv!(recv_buf, comm; source = r, tag = tag))
    end

    # --- Post sends grouped by destination rank (only solid mirror cells) ---
    for r in 0:mpisize-1
        r == myrank && continue
        mp_start = Int(mirror_proc_offsets[r+1]) + 1
        mp_end   = Int(mirror_proc_offsets[r+2])
        mp_start > mp_end && continue
        send_elems = 0
        for k in mp_start:mp_end
            send_elems += mirror_elem_sizes[Int(mirror_proc_mirrors[k]) + 1]
        end
        send_elems == 0 && continue
        send_buf = Vector{Float64}(undef, send_elems)
        off = 0
        for k in mp_start:mp_end
            m_idx = Int(mirror_proc_mirrors[k]) + 1
            sz = mirror_elem_sizes[m_idx]
            sz == 0 && continue
            copyto!(send_buf, off + 1, mirror_data_bufs[m_idx], 1, sz)
            off += sz
        end
        push!(send_bufs, send_buf)
        push!(reqs, MPI.Isend(send_buf, comm; dest = r, tag = tag))
    end
    MPI.Waitall(reqs)

    # --- Scatter received data into the correct (non-contiguous) positions ---
    for (recv_buf, solid_idx, solid_szs) in recv_infos
        off = 0
        for (i, sz) in zip(solid_idx, solid_szs)
            dst = @view ghost_buffer[ghost_offsets[i]+1 : ghost_offsets[i]+sz]
            copyto!(dst, 1, recv_buf, off + 1, sz)
            off += sz
        end
    end
    return nothing
end

"""
$(TYPEDSIGNATURES)
Update `df` in [`GhostVsData`](@ref) for immersed-boundary cells only.

Only ghost cells with `bound_enc < 0` (solid-adjacent, non-InsideSolid) are
communicated; all other ghost cells are skipped, reducing both message sizes
and the number of communicated elements compared to a full [`data_exchange!`](@ref).
"""
function solid_exchange!(p4est::P_pxest_t, ka::KA{DIM,NDF}) where {DIM,NDF}
    MPI.Comm_size(MPI.COMM_WORLD) == 1 && return nothing
    isempty(ka.kinfo.config.IB) && return nothing

    update_solid_mirror_data!(p4est, ka)

    gb_buf     = ka.kdata.ghost.ghost_buffer
    gi         = ka.kdata.ghost.ghost_info
    ghost_wrap = ka.kdata.ghost.ghost_wrap

    # Ghost element sizes: non-zero only for solid-adjacent (bound_enc < 0) ghost cells.
    ghost_szs = [(!isa(gw, GhostInsideSolidData) && gw.bound_enc < 0) ?
                 size_Ghost_Data(gi.ghost_vsnums[i], DIM, NDF) : 0
                 for (i, gw) in enumerate(ghost_wrap)]

    # Mirror element sizes: non-zero only for solid-adjacent mirror cells.
    n_mirrors  = length(gi.mirror_vsnums)
    mirror_szs = Vector{Int}(undef, n_mirrors)
    GC.@preserve p4est ka begin
        ghost_pw = PointerWrapper(ka.kinfo.forest.ghost)
        pp = PointerWrapper(p4est)
        for i in 1:n_mirrors
            pq = pw_mirror_quadrant(pp, ghost_pw, i)
            dp = PointerWrapper(P4estPsData, pq.p.user_data[])
            ps_data = unsafe_pointer_to_objref(pointer(dp.ps_data))
            is_solid = !isa(ps_data, InsideSolidData) && ps_data.bound_enc < 0
            mirror_szs[i] = is_solid ? size_Ghost_Data(gi.mirror_vsnums[i], DIM, NDF) : 0
        end
    end

    _mpi_solid_exchange_varsize!(
        ka.kinfo.forest.ghost,
        gb_buf.ghost_datas, gi.ghost_data_offsets,
        ghost_szs, gb_buf.mirror_data_bufs, mirror_szs, 51)
    return nothing
end

"""
$(TYPEDSIGNATURES)
Non-blocking start of `solid_exchange!`.  Posts `MPI.Irecv!` and `MPI.Isend`
for solid-adjacent ghost cells and returns `recv_bufs` that must be kept alive
until [`solid_exchange_finish!`](@ref) is called.

MPI request handles are stored in `ka.kinfo.status.mpi_reqs`.
"""
function solid_exchange_begin!(p4est::P_pxest_t, ka::KA{DIM,NDF}) where {DIM,NDF}
    reqs = ka.kinfo.status.mpi_reqs
    empty!(reqs)
    recv_bufs = Tuple{Vector{Float64}, Vector{Int}, Vector{Int}}[]

    MPI.Comm_size(MPI.COMM_WORLD) == 1 && return recv_bufs
    isempty(ka.kinfo.config.IB) && return recv_bufs

    update_solid_mirror_data!(p4est, ka)

    gb_buf     = ka.kdata.ghost.ghost_buffer
    gi         = ka.kdata.ghost.ghost_info
    ghost_wrap = ka.kdata.ghost.ghost_wrap

    ghost_szs = [(!isa(gw, GhostInsideSolidData) && gw.bound_enc < 0) ?
                 size_Ghost_Data(gi.ghost_vsnums[i], DIM, NDF) : 0
                 for (i, gw) in enumerate(ghost_wrap)]

    n_mirrors  = length(gi.mirror_vsnums)
    mirror_szs = Vector{Int}(undef, n_mirrors)
    GC.@preserve p4est ka begin
        ghost_pw = PointerWrapper(ka.kinfo.forest.ghost)
        pp = PointerWrapper(p4est)
        for i in 1:n_mirrors
            pq = pw_mirror_quadrant(pp, ghost_pw, i)
            dp = PointerWrapper(P4estPsData, pq.p.user_data[])
            ps_data = unsafe_pointer_to_objref(pointer(dp.ps_data))
            is_solid = !isa(ps_data, InsideSolidData) && ps_data.bound_enc < 0
            mirror_szs[i] = is_solid ? size_Ghost_Data(gi.mirror_vsnums[i], DIM, NDF) : 0
        end
    end

    ghost   = ka.kinfo.forest.ghost
    mpisize, proc_offsets, mirror_proc_offsets, mirror_proc_mirrors =
        _ghost_comm_arrays(ghost)
    myrank = MPI.Comm_rank(MPI.COMM_WORLD)
    comm   = MPI.COMM_WORLD

    # Post receives into per-rank temporary buffers
    for r in 0:mpisize-1
        r == myrank && continue
        g_start = Int(proc_offsets[r+1]) + 1
        g_end   = Int(proc_offsets[r+2])
        g_start > g_end && continue
        solid_idx  = Int[]
        solid_szs  = Int[]
        recv_elems = 0
        for i in g_start:g_end
            sz = ghost_szs[i]
            sz == 0 && continue
            push!(solid_idx, i)
            push!(solid_szs, sz)
            recv_elems += sz
        end
        recv_elems == 0 && continue
        recv_buf = Vector{Float64}(undef, recv_elems)
        push!(recv_bufs, (recv_buf, solid_idx, solid_szs))
        push!(reqs, MPI.Irecv!(recv_buf, comm; source = r, tag = 51))
    end

    # Post sends
    send_bufs = Vector{Vector{Float64}}(undef, 0)
    for r in 0:mpisize-1
        r == myrank && continue
        mp_start = Int(mirror_proc_offsets[r+1]) + 1
        mp_end   = Int(mirror_proc_offsets[r+2])
        mp_start > mp_end && continue
        send_elems = 0
        for k in mp_start:mp_end
            send_elems += mirror_szs[Int(mirror_proc_mirrors[k]) + 1]
        end
        send_elems == 0 && continue
        send_buf = Vector{Float64}(undef, send_elems)
        off = 0
        for k in mp_start:mp_end
            m_idx = Int(mirror_proc_mirrors[k]) + 1
            sz = mirror_szs[m_idx]
            sz == 0 && continue
            copyto!(send_buf, off + 1, gb_buf.mirror_data_bufs[m_idx], 1, sz)
            off += sz
        end
        push!(send_bufs, send_buf)
        push!(reqs, MPI.Isend(send_buf, comm; dest = r, tag = 51))
    end

    return recv_bufs
end

"""
$(TYPEDSIGNATURES)
Complete the asynchronous `solid_exchange!`.  Calls `MPI.Waitall` on the
requests stored in `ka.kinfo.status.mpi_reqs`, then scatters received data
into `ghost_datas`.

`recv_bufs` is the value returned by [`solid_exchange_begin!`](@ref).
"""
function solid_exchange_finish!(ka::KA{DIM,NDF}, recv_bufs) where {DIM,NDF}
    MPI.Comm_size(MPI.COMM_WORLD) == 1 && return nothing
    isempty(ka.kinfo.config.IB) && return nothing

    reqs = ka.kinfo.status.mpi_reqs
    MPI.Waitall(reqs)

    gi = ka.kdata.ghost.ghost_info
    ghost_buffer = ka.kdata.ghost.ghost_buffer.ghost_datas
    ghost_offsets = gi.ghost_data_offsets
    for (recv_buf, solid_idx, solid_szs) in recv_bufs
        off = 0
        for (i, sz) in zip(solid_idx, solid_szs)
            dst = @view ghost_buffer[ghost_offsets[i]+1 : ghost_offsets[i]+sz]
            copyto!(dst, 1, recv_buf, off + 1, sz)
            off += sz
        end
    end
    return nothing
end
