size_Ghost_Data(vs_num,DIM,NDF) = 3 * DIM + 4 + NDF * vs_num
size_Ghost_Slope(vs_num,DIM,NDF) = vs_num * DIM * NDF + DIM*(DIM+2)
size_Ghost_SW(DIM) = DIM * (DIM + 2)
size_Ghost_VS_Structure(vs_num,DIM) = vs_num * (DIM + 2)

"""
$(TYPEDSIGNATURES)
Get the globally largest number of the velocity cells among all local and ghost quadrants.
"""
function get_vs_num(forest::P_pxest_t, ghost::P_pxest_ghost_t)
    get_vs_num(PointerWrapper(forest), PointerWrapper(ghost))
end
function pw_mirror_quadrant(pp::PointerWrapper{p4est_t},gp::PointerWrapper{p4est_ghost_t},i::Integer)
        pm = iPointerWrapper(gp.mirrors, p4est_quadrant_t, i - 1)
        pt = iPointerWrapper(pp.trees, p4est_tree_t, pm.p.piggy3.which_tree[])
        pq = iPointerWrapper(
            pt.quadrants,
            p4est_quadrant_t,
            pm.p.piggy3.local_num[] - pt.quadrants_offset[],
        )
end
function pw_mirror_quadrant(pp::PointerWrapper{p8est_t},gp::PointerWrapper{p8est_ghost_t},i::Integer)
    pm = iPointerWrapper(gp.mirrors, p8est_quadrant_t, i - 1)
    pt = iPointerWrapper(pp.trees, p8est_tree_t, pm.p.piggy3.which_tree[])
    pq = iPointerWrapper(
        pt.quadrants,
        p8est_quadrant_t,
        pm.p.piggy3.local_num[] - pt.quadrants_offset[],
    )
end
function get_vs_num(pp::PW_pxest_t, gp::PW_pxest_ghost_t)
    vs_num = 0
    for i = 1:gp.mirrors.elem_count[]
        pq = pw_mirror_quadrant(pp,gp,i)
        dp = PointerWrapper(P4estPsData, pq.p.user_data[])
        ps_data = unsafe_pointer_to_objref(pointer(dp.ps_data))
        isa(ps_data,InsideSolidData)&&continue
        vs_num = max(vs_num, ps_data.vs_data.vs_num)
    end
    vs_num = MPI.Allreduce(vs_num, max, MPI.COMM_WORLD)
end
function get_mirror_data_inner!(ps_data::PsData{DIM,NDF}, vs_temp::AbstractVector) where{DIM,NDF}
    vs_data = ps_data.vs_data
    vs_num = vs_data.vs_num
    vs_temp[1:NDF*vs_num] .= reshape(vs_data.df, vs_num * NDF)
end
function get_mirror_slope_inner!(ps_data::PsData{DIM,NDF}, vs_temp::AbstractVector) where{DIM,NDF}
    vs_data = ps_data.vs_data
    vs_num = size(vs_data.weight, 1)
    vs_temp[1:(vs_num*NDF*DIM)] .= reshape(vs_data.sdf, vs_num * NDF * DIM)
    vs_temp[vs_num*NDF*DIM+1:vs_num*NDF*DIM+DIM*(DIM+2)] .= reshape(ps_data.sw, DIM*(DIM+2))
end
function get_mirror_slope_inner!(::InsideSolidData, vs_temp::AbstractVector)
    vs_temp[1] = EPS
end
function get_mirror_structure_inner!(ps_data::PsData{DIM}, weight_temp, level_temp, midpoint_temp) where{DIM}
    vs_data = ps_data.vs_data
    vs_num = vs_data.vs_num
    weight_temp[1:vs_num] .= vs_data.weight
    level_temp[1:vs_num] .= vs_data.level
    midpoint_temp[1:(vs_num*DIM)] .= reshape(vs_data.midpoint, vs_num * DIM)
end

function ghost_data_alloc(N::Int, ghost::P_pxest_ghost_t)
    Ptr{Float64}(sc_malloc(P4est.package_id(), sizeof(Float64) * N * PointerWrapper(ghost).ghosts.elem_count[]))
end

# Julia-managed flat ghost receive buffer (replaces sc_malloc'd Ptr version).
function ghost_vec_alloc(total_elems::Int)
    Vector{Float64}(undef, max(total_elems, 1))
end

"""
$(TYPEDSIGNATURES)
Fixed-size ghost exchange (uses `p4est_ghost_exchange_custom`). Kept for the
preliminary size-exchange step in [`exchange_ghost_vsnums`](@ref).
"""
function amr_exchange(
    p4est::Ptr{p4est_t},
    ghost::Ptr{p4est_ghost_t},
    mirror_data_pointers::Array{Ptr{Float64}},
    N::Int,
)
    GC.@preserve mirror_data_pointers ghost N begin
        ghost_datas = ghost_data_alloc(N, ghost)
        GC.@preserve ghost_datas p4est_ghost_exchange_custom(
            p4est,
            ghost,
            sizeof(Float64) * N,
            mirror_data_pointers,
            ghost_datas,
        )
    end
    return ghost_datas
end
function amr_exchange(
    p4est::Ptr{p8est_t},
    ghost::Ptr{p8est_ghost_t},
    mirror_data_pointers::Array{Ptr{Float64}},
    N::Int,
)
    GC.@preserve mirror_data_pointers ghost N begin
        ghost_datas = ghost_data_alloc(N, ghost)
        GC.@preserve ghost_datas p8est_ghost_exchange_custom(
            p4est,
            ghost,
            sizeof(Float64) * N,
            mirror_data_pointers,
            ghost_datas,
        )
    end
    return ghost_datas
end

function ghost_data_ptr(N::Int, ghost_data::Ptr{T}, qid::Integer) where {T}
    return Ptr{Float64}(ghost_data + qid * sizeof(T) * N)
end

# ---------------------------------------------------------------------------
# Variable-size ghost exchange via direct MPI communication
# ---------------------------------------------------------------------------

"""
Extract the p4est ghost communication topology arrays from the ghost pointer.
Returns `(mpisize, proc_offsets, mirror_proc_offsets, mirror_proc_mirrors)`.

Julia indexing convention (1-based):
- Ghosts from rank `r` occupy ghost indices `proc_offsets[r+1]+1 : proc_offsets[r+2]`.
- Mirrors going to rank `r` are at mirror_proc_mirrors indices
  `mirror_proc_offsets[r+1]+1 : mirror_proc_offsets[r+2]`, where each entry is a
  0-based mirror index (add 1 for Julia indexing).
"""
function _ghost_comm_arrays(ghost::P_pxest_ghost_t)
    gp = PointerWrapper(ghost)
    mpisize = Int(gp.mpisize[])
    proc_offsets = Base.unsafe_wrap(
        Vector{Int32}, pointer(gp.proc_offsets), mpisize + 1, own = false)
    mpo = Base.unsafe_wrap(
        Vector{Int32}, pointer(gp.mirror_proc_offsets), mpisize + 1, own = false)
    n_mpm = Int(mpo[end])
    mpm = n_mpm > 0 ?
        Base.unsafe_wrap(Vector{Int32}, pointer(gp.mirror_proc_mirrors), n_mpm, own = false) :
        Int32[]
    return mpisize, proc_offsets, mpo, mpm
end

"""
$(TYPEDSIGNATURES)
Exchange `vs_num` for every quadrant across the ghost layer using a single
fixed-size (`N = 1`) call to `p4est_ghost_exchange_custom`.  Also updates
`kinfo.status.max_vs_num` via `MPI.Allreduce`.

Returns `(ghost_vsnums, mirror_vsnums)`:
- `ghost_vsnums[i]` — `vs_num` of ghost quadrant `i` (0 for `InsideSolidData`).
- `mirror_vsnums[i]` — `vs_num` of mirror quadrant `i` (computed locally).
"""
function exchange_ghost_vsnums(p4est, kinfo::KInfo{DIM,NDF}) where {DIM,NDF}
    GC.@preserve p4est kinfo begin
        gp = PointerWrapper(kinfo.forest.ghost)
        pp = PointerWrapper(p4est)
        n_mirrors = Int(gp.mirrors.elem_count[])
        n_ghosts  = Int(gp.ghosts.elem_count[])

        mirror_vsnums = Vector{Int}(undef, n_mirrors)
        size_ptrs = Array{Ptr{Float64}}(undef, n_mirrors)
        for i in 1:n_mirrors
            pq = pw_mirror_quadrant(pp, gp, i)
            dp = PointerWrapper(P4estPsData, pq.p.user_data[])
            ps_data = unsafe_pointer_to_objref(pointer(dp.ps_data))
            vsn = isa(ps_data, InsideSolidData) ? 0 : ps_data.vs_data.vs_num
            mirror_vsnums[i] = vsn
            p = Ptr{Float64}(sc_malloc(P4est.package_id(), sizeof(Float64)))
            Base.unsafe_store!(p, Float64(vsn))
            size_ptrs[i] = p
        end

        ghost_raw = amr_exchange(p4est, kinfo.forest.ghost, size_ptrs, 1)

        ghost_vsnums = Vector{Int}(undef, n_ghosts)
        for i in 1:n_ghosts
            p = ghost_data_ptr(1, ghost_raw, i - 1)
            ghost_vsnums[i] = Int(round(Base.unsafe_load(p)))
        end

        # local_max = isempty(mirror_vsnums) ? 0 : maximum(mirror_vsnums)
        # ghost_max = isempty(ghost_vsnums)  ? 0 : maximum(ghost_vsnums)
        # kinfo.status.max_vs_num = Int(MPI.Allreduce(max(local_max, ghost_max), max, MPI.COMM_WORLD))

        id = P4est.package_id()
        for p in size_ptrs; sc_free(id, p); end
        sc_free(id, ghost_raw)
    end
    return ghost_vsnums, mirror_vsnums
end

"""
Core variable-size MPI exchange. Chunks are packed in the native p4est ghost
topology order: receiver ghosts use `proc_offsets`, sender mirrors use
`mirror_proc_mirrors`. Zero-sized entries are skipped in that same order, which
lets the same routine handle level-filtered exchanges while keeping the fixed
ghost offsets.
"""
function _mpi_exchange_varsize!(
    ghost::P_pxest_ghost_t,
    ghost_buffer::Vector{Float64},
    ghost_offsets::Vector{Int},
    ghost_elem_sizes::Vector{Int},
    mirror_data_bufs::Vector{Vector{Float64}},
    mirror_elem_sizes::Vector{Int},
    tag::Int,
)
    mpisize, proc_offsets, mirror_proc_offsets, mirror_proc_mirrors =
        _ghost_comm_arrays(ghost)
    myrank = MPI.Comm_rank(MPI.COMM_WORLD)
    comm = MPI.COMM_WORLD

    reqs       = Vector{MPI.Request}(undef, 0)
    send_bufs  = Vector{Vector{Float64}}(undef, 0)
    recv_infos = Vector{Tuple{Any, Vector{Int}, Bool}}(undef, 0)

    for r in 0:mpisize-1
        r == myrank && continue
        ghost_ids = _active_ghost_ids_native(ghost_elem_sizes, proc_offsets, r)
        isempty(ghost_ids) && continue
        recv_elems = 0
        for i in ghost_ids
            recv_elems += ghost_elem_sizes[i]
        end
        recv_elems == 0 && continue
        if _is_compact_recv_layout(ghost_ids, ghost_offsets, ghost_elem_sizes)
            first_offset = ghost_offsets[first(ghost_ids)]
            recv_buf = @view ghost_buffer[first_offset+1:first_offset+recv_elems]
            push!(recv_infos, (recv_buf, Int[], false))
        else
            recv_buf = Vector{Float64}(undef, recv_elems)
            push!(recv_infos, (recv_buf, ghost_ids, true))
        end
        push!(reqs, MPI.Irecv!(recv_buf, comm; source = r, tag = tag))
    end
    n_recvs = length(recv_infos)

    for r in 0:mpisize-1
        r == myrank && continue
        mp_start = Int(mirror_proc_offsets[r+1]) + 1
        mp_end   = Int(mirror_proc_offsets[r+2])
        mp_start > mp_end && continue
        send_elems = 0
        for k in mp_start:mp_end
            m_idx = Int(mirror_proc_mirrors[k]) + 1
            send_elems += mirror_elem_sizes[m_idx]
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

    isempty(reqs) && return nothing
    statuses = GC.@preserve send_bufs recv_infos MPI.Waitall(reqs, MPI.Status)

    @inbounds for recv_id in 1:n_recvs
        recv_buf, ghost_ids, needs_scatter = recv_infos[recv_id]
        recv_count = MPI.Get_count(statuses[recv_id], Float64)
        recv_count == length(recv_buf) || error(
            "Direct MPI ghost exchange size mismatch on rank ", myrank,
            ": expected ", length(recv_buf), " Float64 values, received ", recv_count)
        needs_scatter || continue
        off = 0
        for i in ghost_ids
            sz = ghost_elem_sizes[i]
            copyto!(ghost_buffer, ghost_offsets[i] + 1, recv_buf, off + 1, sz)
            off += sz
        end
    end
    return nothing
end

"""
$(TYPEDSIGNATURES)
Sparse variant of [`_mpi_exchange_varsize!`](@ref) for level-filtered
exchanges where some entries of `ghost_elem_sizes` / `mirror_elem_sizes`
are zero.
"""
function _mpi_exchange_varsize_sparse!(
    ghost::P_pxest_ghost_t,
    ghost_buffer::Vector{Float64},
    ghost_offsets::Vector{Int},
    ghost_elem_sizes::Vector{Int},
    mirror_data_bufs::Vector{Vector{Float64}},
    mirror_elem_sizes::Vector{Int},
    tag::Int,
)
    _mpi_exchange_varsize!(ghost, ghost_buffer, ghost_offsets, ghost_elem_sizes,
                           mirror_data_bufs, mirror_elem_sizes, tag)
end

function _active_ghost_ids_native(ghost_elem_sizes::Vector{Int}, proc_offsets, r::Integer)
    g_start = Int(proc_offsets[r+1]) + 1
    g_end   = Int(proc_offsets[r+2])
    ids = Int[]
    g_start > g_end && return ids
    @inbounds for i in g_start:g_end
        ghost_elem_sizes[i] == 0 && continue
        push!(ids, i)
    end
    return ids
end

function _is_compact_recv_layout(ids::Vector{Int}, ghost_offsets::Vector{Int},
                                 ghost_elem_sizes::Vector{Int})
    isempty(ids) && return true
    expected = ghost_offsets[first(ids)]
    @inbounds for i in ids
        ghost_offsets[i] == expected || return false
        expected += ghost_elem_sizes[i]
    end
    return true
end

function _native_exchange_varsize(
    p4est::P_pxest_t,
    ghost::P_pxest_ghost_t,
    mirror_data_bufs::Vector{Vector{Float64}},
    mirror_elem_sizes::Vector{Int},
    ghost_elem_sizes::Vector{Int},
)
    max_mirror = isempty(mirror_elem_sizes) ? 0 : maximum(mirror_elem_sizes)
    max_ghost = isempty(ghost_elem_sizes) ? 0 : maximum(ghost_elem_sizes)
    local_N = max(1, max_mirror, max_ghost)
    N = MPI.Allreduce(local_N, max, MPI.COMM_WORLD)

    send_bufs = [zeros(Float64, N) for _ in eachindex(mirror_data_bufs)]
    ptrs = Vector{Ptr{Float64}}(undef, length(send_bufs))
    @inbounds for i in eachindex(send_bufs)
        sz = mirror_elem_sizes[i]
        sz > 0 && copyto!(send_bufs[i], 1, mirror_data_bufs[i], 1, sz)
        ptrs[i] = pointer(send_bufs[i])
    end

    ghost_raw = GC.@preserve send_bufs ptrs amr_exchange(p4est, ghost, ptrs, N)

    n_ghosts = length(ghost_elem_sizes)
    ghost_offsets = Vector{Int}(undef, n_ghosts)
    total_elems = 0
    @inbounds for i in 1:n_ghosts
        ghost_offsets[i] = total_elems
        total_elems += ghost_elem_sizes[i]
    end
    ghost_buffer = ghost_vec_alloc(total_elems)
    @inbounds for i in 1:n_ghosts
        sz = ghost_elem_sizes[i]
        sz == 0 && continue
        src = ghost_data_ptr(N, ghost_raw, i - 1)
        src_vec = Base.unsafe_wrap(Vector{Float64}, src, N, own = false)
        copyto!(ghost_buffer, ghost_offsets[i] + 1, src_vec, 1, sz)
    end
    sc_free(P4est.package_id(), ghost_raw)
    return ghost_buffer, ghost_offsets
end

function _native_exchange_varsize!(
    p4est::P_pxest_t,
    ghost::P_pxest_ghost_t,
    ghost_buffer::Vector{Float64},
    ghost_offsets::Vector{Int},
    ghost_elem_sizes::Vector{Int},
    mirror_data_bufs::Vector{Vector{Float64}},
    mirror_elem_sizes::Vector{Int},
)
    max_mirror = isempty(mirror_elem_sizes) ? 0 : maximum(mirror_elem_sizes)
    max_ghost = isempty(ghost_elem_sizes) ? 0 : maximum(ghost_elem_sizes)
    local_N = max(1, max_mirror, max_ghost)
    N = MPI.Allreduce(local_N, max, MPI.COMM_WORLD)

    send_bufs = [zeros(Float64, N) for _ in eachindex(mirror_data_bufs)]
    ptrs = Vector{Ptr{Float64}}(undef, length(send_bufs))
    @inbounds for i in eachindex(send_bufs)
        sz = mirror_elem_sizes[i]
        sz > 0 && copyto!(send_bufs[i], 1, mirror_data_bufs[i], 1, sz)
        ptrs[i] = pointer(send_bufs[i])
    end

    ghost_raw = GC.@preserve send_bufs ptrs amr_exchange(p4est, ghost, ptrs, N)
    @inbounds for i in eachindex(ghost_elem_sizes)
        sz = ghost_elem_sizes[i]
        sz == 0 && continue
        src = ghost_data_ptr(N, ghost_raw, i - 1)
        src_vec = Base.unsafe_wrap(Vector{Float64}, src, N, own = false)
        copyto!(ghost_buffer, ghost_offsets[i] + 1, src_vec, 1, sz)
    end
    sc_free(P4est.package_id(), ghost_raw)
    return nothing
end

"""
$(TYPEDSIGNATURES)
Allocate a variable-size ghost receive buffer, perform the MPI exchange, and
return `(ghost_buffer, ghost_offsets)`.

- `ghost_buffer`: Julia-managed `Vector{Float64}` holding all ghost data concatenated.
- `ghost_offsets[i]`: Float64-element offset of ghost `i`'s data within `ghost_buffer`.
"""
function amr_exchange_varsize(
    p4est::P_pxest_t,
    ghost::P_pxest_ghost_t,
    mirror_data_bufs::Vector{Vector{Float64}},
    mirror_elem_sizes::Vector{Int},
    ghost_elem_sizes::Vector{Int},
    tag::Int,
)
    ghost_buffer, ghost_offsets =
        amr_exchange_varsize(ghost, mirror_data_bufs, mirror_elem_sizes, ghost_elem_sizes, tag)
    return ghost_buffer, ghost_offsets
end

function amr_exchange_varsize(
    ghost::P_pxest_ghost_t,
    mirror_data_bufs::Vector{Vector{Float64}},
    mirror_elem_sizes::Vector{Int},
    ghost_elem_sizes::Vector{Int},
    tag::Int,
)
    n_ghosts = length(ghost_elem_sizes)
    ghost_offsets = Vector{Int}(undef, n_ghosts)
    total_elems = 0
    for i in 1:n_ghosts
        ghost_offsets[i] = total_elems
        total_elems += ghost_elem_sizes[i]
    end
    ghost_buffer = ghost_vec_alloc(total_elems)
    if MPI.Comm_size(MPI.COMM_WORLD) > 1
        _mpi_exchange_varsize!(ghost, ghost_buffer, ghost_offsets, ghost_elem_sizes,
                               mirror_data_bufs, mirror_elem_sizes, tag)
    end
    return ghost_buffer, ghost_offsets
end

"""
$(TYPEDSIGNATURES)
In-place variable-size ghost exchange: re-uses the pre-allocated `ghost_buffer`
(from a previous [`amr_exchange_varsize`](@ref)) and re-exchanges new mirror data
without reallocating.
"""
function amr_exchange_varsize!(
    p4est::P_pxest_t,
    ghost::P_pxest_ghost_t,
    ghost_buffer::Vector{Float64},
    ghost_offsets::Vector{Int},
    ghost_elem_sizes::Vector{Int},
    mirror_data_bufs::Vector{Vector{Float64}},
    mirror_elem_sizes::Vector{Int},
    tag::Int,
)
    amr_exchange_varsize!(
        ghost, ghost_buffer, ghost_offsets, ghost_elem_sizes,
        mirror_data_bufs, mirror_elem_sizes, tag)
    return nothing
end

function amr_exchange_varsize!(
    ghost::P_pxest_ghost_t,
    ghost_buffer::Vector{Float64},
    ghost_offsets::Vector{Int},
    ghost_elem_sizes::Vector{Int},
    mirror_data_bufs::Vector{Vector{Float64}},
    mirror_elem_sizes::Vector{Int},
    tag::Int,
)
    MPI.Comm_size(MPI.COMM_WORLD) == 1 && return nothing
    _mpi_exchange_varsize!(ghost, ghost_buffer, ghost_offsets, ghost_elem_sizes,
                           mirror_data_bufs, mirror_elem_sizes, tag)
    return nothing
end

"""
$(TYPEDSIGNATURES)
Sparse counterpart of [`amr_exchange_varsize!`](@ref). Use this when some
entries of `mirror_elem_sizes` / `ghost_elem_sizes` are zero (e.g.
level-filtered exchanges). The common direct-MPI core packs only positive-size
chunks and scatters them back to the fixed `ghost_offsets`.
"""
function amr_exchange_varsize_sparse!(
    p4est::P_pxest_t,
    ghost::P_pxest_ghost_t,
    ghost_buffer::Vector{Float64},
    ghost_offsets::Vector{Int},
    ghost_elem_sizes::Vector{Int},
    mirror_data_bufs::Vector{Vector{Float64}},
    mirror_elem_sizes::Vector{Int},
    tag::Int,
)
    amr_exchange_varsize_sparse!(
        ghost, ghost_buffer, ghost_offsets, ghost_elem_sizes,
        mirror_data_bufs, mirror_elem_sizes, tag)
    return nothing
end

function amr_exchange_varsize_sparse!(
    ghost::P_pxest_ghost_t,
    ghost_buffer::Vector{Float64},
    ghost_offsets::Vector{Int},
    ghost_elem_sizes::Vector{Int},
    mirror_data_bufs::Vector{Vector{Float64}},
    mirror_elem_sizes::Vector{Int},
    tag::Int,
)
    MPI.Comm_size(MPI.COMM_WORLD) == 1 && return nothing
    _mpi_exchange_varsize_sparse!(ghost, ghost_buffer, ghost_offsets, ghost_elem_sizes,
                                  mirror_data_bufs, mirror_elem_sizes, tag)
    return nothing
end

# ---------------------------------------------------------------------------
# Mirror buffer packing (exact-size allocation per quadrant)
# ---------------------------------------------------------------------------

function get_mirror_data(p4est, kinfo::KInfo{DIM,NDF}, mirror_vsnums::Vector{Int}) where{DIM,NDF}
    GC.@preserve p4est kinfo begin
        gp = PointerWrapper(kinfo.forest.ghost)
        pp = PointerWrapper(p4est)
        n_mirrors = Int(gp.mirrors.elem_count[])
        mirror_data_bufs = Vector{Vector{Float64}}(undef, n_mirrors)
        for i = 1:n_mirrors
            pq = pw_mirror_quadrant(pp,gp,i)
            dp = PointerWrapper(P4estPsData, pq.p.user_data[])
            ps_data = unsafe_pointer_to_objref(pointer(dp.ps_data))
            vs_num = mirror_vsnums[i]
            sz = size_Ghost_Data(vs_num, DIM, NDF)
            ap = Vector{Float64}(undef, sz)
            offset = 0
            if isa(ps_data,InsideSolidData)
                ap[1:(offset+=DIM)] .= ps_data.ds
                ap[offset+1:(offset+=DIM)] .= ps_data.midpoint
                ap[offset+1:(offset+=DIM+2)] .= Inf
                ap[offset+=1] = 0
                ap[offset+=1] = ps_data.bound_enc
            else
                ap[1:(offset+=DIM)] .= ps_data.ds
                ap[offset+1:(offset+=DIM)] .= ps_data.midpoint
                ap[offset+1:(offset+=DIM+2)] .= ps_data.w
                ap[offset+=1] = ps_data.vs_data.vs_num
                ap[offset+=1] = ps_data.bound_enc
                vs_temp = @view(ap[offset+1:(offset+=NDF*vs_num)])
                get_mirror_data_inner!(ps_data, vs_temp)
            end
            mirror_data_bufs[i] = ap
        end
    end
    return mirror_data_bufs
end

function get_mirror_slope(p4est, kinfo::KInfo{DIM,NDF}, mirror_vsnums::Vector{Int}) where{DIM,NDF}
    GC.@preserve p4est kinfo begin
        gp = PointerWrapper(kinfo.forest.ghost)
        pp = PointerWrapper(p4est)
        n_mirrors = Int(gp.mirrors.elem_count[])
        mirror_slope_bufs = Vector{Vector{Float64}}(undef, n_mirrors)
        for i = 1:n_mirrors
            pq = pw_mirror_quadrant(pp,gp,i)
            dp = PointerWrapper(P4estPsData, pq.p.user_data[])
            ps_data = unsafe_pointer_to_objref(pointer(dp.ps_data))
            vs_num = mirror_vsnums[i]
            sz = size_Ghost_Slope(vs_num, DIM, NDF)
            # For InsideSolidData (vs_num=0), size_Ghost_Slope = DIM*(DIM+2) >= 1.
            ap = Vector{Float64}(undef, max(sz, 1))
            get_mirror_slope_inner!(ps_data, ap)
            mirror_slope_bufs[i] = ap
        end
    end
    return mirror_slope_bufs
end

function get_mirror_structure(p4est, kinfo::KInfo{DIM,NDF}, mirror_vsnums::Vector{Int}) where{DIM,NDF}
    GC.@preserve p4est kinfo begin
        gp = PointerWrapper(kinfo.forest.ghost)
        pp = PointerWrapper(p4est)
        n_mirrors = Int(gp.mirrors.elem_count[])
        mirror_structure_bufs = Vector{Vector{Float64}}(undef, n_mirrors)
        for i = 1:n_mirrors
            pq = pw_mirror_quadrant(pp,gp,i)
            dp = PointerWrapper(P4estPsData, pq.p.user_data[])
            ps_data = unsafe_pointer_to_objref(pointer(dp.ps_data))
            if isa(ps_data, InsideSolidData)
                # size_Ghost_VS_Structure(0,DIM)=0; 1-element placeholder.
                # Zero elements are exchanged for InsideSolidData (ghost_vsnums=0).
                mirror_structure_bufs[i] = Vector{Float64}(undef, 1)
                continue
            end
            vs_num = mirror_vsnums[i]
            sz = size_Ghost_VS_Structure(vs_num, DIM)
            ap = Vector{Float64}(undef, max(sz, 1))
            weight_temp   = @view(ap[1:vs_num])
            level_temp    = @view(ap[vs_num+1:2*vs_num])
            midpoint_temp = @view(ap[2*vs_num+1:vs_num*(DIM+2)])
            get_mirror_structure_inner!(ps_data, weight_temp, level_temp, midpoint_temp)
            mirror_structure_bufs[i] = ap
        end
        return mirror_structure_bufs
    end
end

# ---------------------------------------------------------------------------
# Ghost pointer and wrap initialisation
# ---------------------------------------------------------------------------

function initialize_ghost_pool(p4est, kinfo::KInfo{DIM,NDF}) where{DIM,NDF}
    ghost_vsnums, mirror_vsnums = exchange_ghost_vsnums(p4est, kinfo)
    ghost_levels, mirror_levels = exchange_ghost_levels(p4est, kinfo)

    mirror_data_szs  = [size_Ghost_Data(v, DIM, NDF) for v in mirror_vsnums]
    ghost_data_szs   = [size_Ghost_Data(v, DIM, NDF) for v in ghost_vsnums]
    mirror_data_bufs = get_mirror_data(p4est, kinfo, mirror_vsnums)
    ghost_datas, ghost_data_offsets = amr_exchange_varsize(
        p4est, kinfo.forest.ghost, mirror_data_bufs, mirror_data_szs, ghost_data_szs, 51)

    mirror_slope_szs  = [size_Ghost_Slope(v, DIM, NDF) for v in mirror_vsnums]
    ghost_slope_szs   = [size_Ghost_Slope(v, DIM, NDF) for v in ghost_vsnums]
    mirror_slope_bufs = get_mirror_slope(p4est, kinfo, mirror_vsnums)
    ghost_slopes, ghost_slope_offsets = amr_exchange_varsize(
        p4est, kinfo.forest.ghost, mirror_slope_bufs, mirror_slope_szs, ghost_slope_szs, 52)

    mirror_structure_bufs = get_mirror_structure(p4est, kinfo, mirror_vsnums)
    mirror_struct_szs = [size_Ghost_VS_Structure(v, DIM) for v in mirror_vsnums]
    ghost_struct_szs  = [size_Ghost_VS_Structure(v, DIM) for v in ghost_vsnums]
    ghost_structures, ghost_structure_offsets = amr_exchange_varsize(
        p4est, kinfo.forest.ghost, mirror_structure_bufs, mirror_struct_szs, ghost_struct_szs, 53)

    gb = GhostBuffer(
        ghost_datas, ghost_slopes, ghost_structures,
        mirror_data_bufs, mirror_slope_bufs, mirror_structure_bufs,
    )
    gi = GhostInfo(
        ghost_vsnums, mirror_vsnums,
        ghost_data_offsets, ghost_slope_offsets, ghost_structure_offsets,
        ghost_data_szs, ghost_slope_szs, mirror_data_szs, mirror_slope_szs,
        ghost_levels, mirror_levels,
    )
    return gb, gi
end

"""
$(TYPEDSIGNATURES)
Exchange the physical-space refinement level of every quadrant across the
ghost layer using a fixed-size (`N = 1`) call to `p4est_ghost_exchange_custom`.
Returns `(ghost_levels::Vector{Int8}, mirror_levels::Vector{Int8})`.
"""
function exchange_ghost_levels(p4est, kinfo::KInfo)
    GC.@preserve p4est kinfo begin
        gp = PointerWrapper(kinfo.forest.ghost)
        pp = PointerWrapper(p4est)
        n_mirrors = Int(gp.mirrors.elem_count[])
        n_ghosts  = Int(gp.ghosts.elem_count[])

        mirror_levels = Vector{Int8}(undef, n_mirrors)
        lvl_ptrs = Array{Ptr{Float64}}(undef, n_mirrors)
        for i in 1:n_mirrors
            pq = pw_mirror_quadrant(pp, gp, i)
            lvl = Int8(pq.level[])
            mirror_levels[i] = lvl
            p = Ptr{Float64}(sc_malloc(P4est.package_id(), sizeof(Float64)))
            Base.unsafe_store!(p, Float64(lvl))
            lvl_ptrs[i] = p
        end

        ghost_raw = amr_exchange(p4est, kinfo.forest.ghost, lvl_ptrs, 1)

        ghost_levels = Vector{Int8}(undef, n_ghosts)
        for i in 1:n_ghosts
            p = ghost_data_ptr(1, ghost_raw, i - 1)
            ghost_levels[i] = Int8(round(Base.unsafe_load(p)))
        end

        id = P4est.package_id()
        for p in lvl_ptrs; sc_free(id, p); end
        sc_free(id, ghost_raw)
    end
    return ghost_levels, mirror_levels
end

function initialize_ghost_wrap(kinfo::KInfo{DIM,NDF},
                               ghost_buffer::GhostBuffer,
                               ghost_info::GhostInfo) where{DIM,NDF}
    gp = PointerWrapper(kinfo.forest.ghost)
    n_ghosts = Int(gp.ghosts.elem_count[])
    ghost_wrap = Array{AbstractGhostPsData}(undef, n_ghosts)

    ghost_vsnums            = ghost_info.ghost_vsnums
    ghost_data_offsets      = ghost_info.ghost_data_offsets
    ghost_slope_offsets     = ghost_info.ghost_slope_offsets
    ghost_structure_offsets = ghost_info.ghost_structure_offsets

    # Raw pointers into Julia-managed buffers.  The buffers live in ghost_buffer
    # (held by Ghost → ka.kdata.ghost) so they outlive all the wrapped arrays.
    p_datas   = pointer(ghost_buffer.ghost_datas)
    p_slopes  = pointer(ghost_buffer.ghost_slopes)
    p_structs = pointer(ghost_buffer.ghost_structures)

    for i in eachindex(ghost_wrap)
        pq = iPointerWrapper(gp.ghosts, p4est_quadrant_t, i - 1)
        owner_rank = pq.p.piggy1.owner_rank[]
        which_tree = pq.p.piggy3.which_tree[]
        local_num  = pq.p.piggy3.local_num[]
        quadid = which_tree * 2^(DIM * kinfo.config.solver.AMR_PS_MAXLEVEL) + local_num

        vs_num = ghost_vsnums[i]

        # --- data buffer (ds, midpoint, w, vs_num, bound_enc, df) ---
        p_data = p_datas + ghost_data_offsets[i] * sizeof(Float64)
        offset = 0
        ds       = Base.unsafe_wrap(Vector{Float64}, p_data + offset * sizeof(Float64), DIM)
        offset  += DIM
        midpoint = Base.unsafe_wrap(Vector{Float64}, p_data + offset * sizeof(Float64), DIM)
        offset  += DIM
        w        = Base.unsafe_wrap(Vector{Float64}, p_data + offset * sizeof(Float64), DIM + 2)
        offset  += DIM + 2
        bound_enc = Int(Base.unsafe_load(p_data + (offset + 1) * sizeof(Float64)))

        if w[1] == Inf
            ghost_wrap[i] = GhostInsideSolidData{DIM,NDF}(bound_enc, midpoint, ds)
            continue
        end

        offset += 2
        df = Base.unsafe_wrap(Matrix{Float64}, p_data + offset * sizeof(Float64), (vs_num, NDF))

        # --- slope buffer (sdf, sw) ---
        p_slope = p_slopes + ghost_slope_offsets[i] * sizeof(Float64)
        sdf = Base.unsafe_wrap(Array{Float64},  p_slope, (vs_num, NDF, DIM))
        sw  = Base.unsafe_wrap(Matrix{Float64}, p_slope + vs_num * NDF * DIM * sizeof(Float64),
                               (DIM + 2, DIM))

        # --- structure buffer (weight, level, midpoint_vs) ---
        p_struct    = p_structs + ghost_structure_offsets[i] * sizeof(Float64)
        weight      = Base.unsafe_wrap(Vector{Float64}, p_struct, vs_num)
        level       = Base.unsafe_wrap(Vector{Float64}, p_struct + vs_num * sizeof(Float64), vs_num)
        midpoint_vs = Base.unsafe_wrap(Matrix{Float64},
                                       p_struct + 2 * vs_num * sizeof(Float64), (vs_num, DIM))

        vs_data = GhostVsData{DIM,NDF}(vs_num, Int8.(round.(level)), weight, midpoint_vs, df, sdf)
        ghost_wrap[i] = GhostPsData{DIM,NDF}(owner_rank, quadid, bound_enc, ds, midpoint, w, sw, vs_data)
    end
    return ghost_wrap
end

# ---------------------------------------------------------------------------
# Mirror buffer update (in-place, called each time-step before exchange)
# ---------------------------------------------------------------------------

function update_mirror_data!(p4est, ka::KA{DIM,NDF}) where{DIM,NDF}
    GC.@preserve p4est ka begin
        mirror_data_bufs = ka.kdata.ghost.ghost_buffer.mirror_data_bufs
        mirror_vsnums    = ka.kdata.ghost.ghost_info.mirror_vsnums
        gp = PointerWrapper(ka.kinfo.forest.ghost)
        pp = PointerWrapper(p4est)
        for i in eachindex(mirror_data_bufs)
            pq = pw_mirror_quadrant(pp,gp,i)
            dp = PointerWrapper(P4estPsData, pq.p.user_data[])
            ps_data = unsafe_pointer_to_objref(pointer(dp.ps_data))
            isa(ps_data,InsideSolidData) && continue
            ap = mirror_data_bufs[i]
            ap[DIM*2+1:DIM*3+2] .= ps_data.w
            vs_temp = @view(ap[3*DIM+4+1:end])
            get_mirror_data_inner!(ps_data, vs_temp)
        end
    end
end

function update_solid_mirror_data!(p4est, ka::KA{DIM,NDF}) where{DIM,NDF}
    GC.@preserve p4est ka begin
        mirror_data_bufs = ka.kdata.ghost.ghost_buffer.mirror_data_bufs
        mirror_vsnums    = ka.kdata.ghost.ghost_info.mirror_vsnums
        gp = PointerWrapper(ka.kinfo.forest.ghost)
        pp = PointerWrapper(p4est)
        for i in eachindex(mirror_data_bufs)
            pq = pw_mirror_quadrant(pp,gp,i)
            dp = PointerWrapper(P4estPsData, pq.p.user_data[])
            ps_data = unsafe_pointer_to_objref(pointer(dp.ps_data))
            (isa(ps_data,InsideSolidData) || ps_data.bound_enc >= 0) && continue
            ap = mirror_data_bufs[i]
            ap[DIM*2+1:DIM*3+2] .= ps_data.w
            vs_temp = @view(ap[3*DIM+4+1:end])
            get_mirror_data_inner!(ps_data, vs_temp)
        end
    end
end

function update_mirror_slope!(p4est, ka::KA{DIM,NDF}) where{DIM,NDF}
    GC.@preserve p4est ka begin
        mirror_slope_bufs = ka.kdata.ghost.ghost_buffer.mirror_slope_bufs
        gp = PointerWrapper(ka.kinfo.forest.ghost)
        pp = PointerWrapper(p4est)
        for i in eachindex(mirror_slope_bufs)
            pq = pw_mirror_quadrant(pp,gp,i)
            dp = PointerWrapper(P4estPsData, pq.p.user_data[])
            ps_data = unsafe_pointer_to_objref(pointer(dp.ps_data))
            isa(ps_data,InsideSolidData) && continue
            get_mirror_slope_inner!(ps_data, mirror_slope_bufs[i])
        end
    end
end

function get_mirror_sw(p4est, ka::KA{DIM,NDF}) where{DIM,NDF}
    sw_sz = size_Ghost_SW(DIM)
    GC.@preserve p4est ka begin
        gp = PointerWrapper(ka.kinfo.forest.ghost)
        pp = PointerWrapper(p4est)
        n_mirrors = Int(gp.mirrors.elem_count[])
        mirror_sw_bufs = Vector{Vector{Float64}}(undef, n_mirrors)
        for i in 1:n_mirrors
            pq = pw_mirror_quadrant(pp, gp, i)
            dp = PointerWrapper(P4estPsData, pq.p.user_data[])
            ps_data = unsafe_pointer_to_objref(pointer(dp.ps_data))
            buf = Vector{Float64}(undef, sw_sz)
            if isa(ps_data, InsideSolidData)
                fill!(buf, 0.0)
            else
                buf .= reshape(ps_data.sw, sw_sz)
            end
            mirror_sw_bufs[i] = buf
        end
    end
    return mirror_sw_bufs
end

# ---------------------------------------------------------------------------
# High-level exchange functions (public API, called from Solver)
# ---------------------------------------------------------------------------

"""
$(TYPEDSIGNATURES)
Update `df` in [`GhostVsData`](@ref).
"""
function data_exchange!(p4est::P_pxest_t, ka::KA{DIM,NDF}) where{DIM,NDF}
    MPI.Comm_size(MPI.COMM_WORLD) == 1 && return nothing
    update_mirror_data!(p4est, ka)
    gb_buf = ka.kdata.ghost.ghost_buffer
    gi     = ka.kdata.ghost.ghost_info
    amr_exchange_varsize!(
        p4est, ka.kinfo.forest.ghost,
        gb_buf.ghost_datas, gi.ghost_data_offsets, gi.ghost_data_szs,
        gb_buf.mirror_data_bufs, gi.mirror_data_szs, 51)
end

"""
$(TYPEDSIGNATURES)
Update `sdf` in [`GhostVsData`](@ref).
"""
function slope_exchange!(p4est::P_pxest_t, ka::KA{DIM,NDF}) where{DIM,NDF}
    MPI.Comm_size(MPI.COMM_WORLD) == 1 && return nothing
    update_mirror_slope!(p4est, ka)
    gb_buf = ka.kdata.ghost.ghost_buffer
    gi     = ka.kdata.ghost.ghost_info
    _native_exchange_varsize!(
        p4est, ka.kinfo.forest.ghost,
        gb_buf.ghost_slopes, gi.ghost_slope_offsets, gi.ghost_slope_szs,
        gb_buf.mirror_slope_bufs, gi.mirror_slope_szs)
end

function sw_exchange!(p4est::P_pxest_t, ka::KA{DIM,NDF}) where{DIM,NDF}
    MPI.Comm_size(MPI.COMM_WORLD) == 1 && return nothing
    gb_buf = ka.kdata.ghost.ghost_buffer
    gi     = ka.kdata.ghost.ghost_info
    sw_sz = size_Ghost_SW(DIM)
    mirror_sw_bufs = get_mirror_sw(p4est, ka)
    mirror_sw_szs = fill(sw_sz, length(gi.mirror_slope_szs))
    ghost_sw_szs = fill(sw_sz, length(gi.ghost_slope_szs))
    ghost_sw_offsets = Vector{Int}(undef, length(gi.ghost_slope_offsets))
    @inbounds for i in eachindex(ghost_sw_offsets)
        ghost_sw_offsets[i] = gi.ghost_slope_offsets[i] + gi.ghost_vsnums[i] * NDF * DIM
    end
    amr_exchange_varsize!(
        p4est, ka.kinfo.forest.ghost,
        gb_buf.ghost_slopes, ghost_sw_offsets, ghost_sw_szs,
        mirror_sw_bufs, mirror_sw_szs, 56)
    return nothing
end

"""
$(TYPEDSIGNATURES)
Level-filtered variant of [`slope_exchange!`](@ref): only mirrors whose
physical-space refinement level equals `L` participate in the send, and only
ghosts with `ghost_levels == L` are written into the receive buffer. Other
mirror/ghost slots are temporarily zero-sized for this exchange; the buffer
contents at those slots are left untouched. Used by [`slope!`](@ref) to push
the freshly transverse-corrected sw/sdf of one refinement level to ranks
where that level appears as a coarse neighbour of a finer cell.
"""
function slope_exchange_level!(p4est::P_pxest_t, ka::KA{DIM,NDF}, L::Integer) where{DIM,NDF}
    MPI.Comm_size(MPI.COMM_WORLD) == 1 && return nothing
    gb_buf = ka.kdata.ghost.ghost_buffer
    gi     = ka.kdata.ghost.ghost_info

    GC.@preserve p4est ka begin
        mirror_slope_bufs = gb_buf.mirror_slope_bufs
        gp = PointerWrapper(ka.kinfo.forest.ghost)
        pp = PointerWrapper(p4est)
        @inbounds for i in eachindex(mirror_slope_bufs)
            gi.mirror_levels[i] == L || continue
            pq = pw_mirror_quadrant(pp, gp, i)
            dp = PointerWrapper(P4estPsData, pq.p.user_data[])
            ps_data = unsafe_pointer_to_objref(pointer(dp.ps_data))
            isa(ps_data, InsideSolidData) && continue
            get_mirror_slope_inner!(ps_data, mirror_slope_bufs[i])
        end
    end

    n_mirrors = length(gi.mirror_slope_szs)
    n_ghosts  = length(gi.ghost_slope_szs)
    mirror_szs = Vector{Int}(undef, n_mirrors)
    ghost_szs  = Vector{Int}(undef, n_ghosts)
    @inbounds for i in 1:n_mirrors
        mirror_szs[i] = (gi.mirror_levels[i] == L) ? gi.mirror_slope_szs[i] : 0
    end
    @inbounds for i in 1:n_ghosts
        ghost_szs[i] = (gi.ghost_levels[i] == L) ? gi.ghost_slope_szs[i] : 0
    end
    amr_exchange_varsize_sparse!(
        p4est, ka.kinfo.forest.ghost,
        gb_buf.ghost_slopes, gi.ghost_slope_offsets, ghost_szs,
        gb_buf.mirror_slope_bufs, mirror_szs, 55)
end

"""
$(TYPEDSIGNATURES)
Exchange a single Bool-as-Float per mirror cell: whether
`maximum(ps_data.lohner) > ADAPT_COEFFI_PS`. Returns a map keyed by ghost-object
identity, with `true` for ghosts whose owner-side cell is above threshold.
Caller must have already filled `ps_data.lohner` on all local cells.
`GhostInsideSolidData` ghosts are left out of the map.
"""
function lohner_flag_exchange!(p4est::P_pxest_t, ka::KA{DIM,NDF}) where{DIM,NDF}
    ghost_flags = Dict{UInt64, Bool}()
    MPI.Comm_size(MPI.COMM_WORLD) == 1 && return ghost_flags
    ghost     = ka.kinfo.forest.ghost
    threshold = ka.kinfo.config.solver.ADAPT_COEFFI_PS

    gp = PointerWrapper(ghost)
    pp = PointerWrapper(p4est)
    n_mirrors = Int(gp.mirrors.elem_count[])
    n_ghosts  = Int(gp.ghosts.elem_count[])

    mirror_bufs = Vector{Vector{Float64}}(undef, n_mirrors)
    GC.@preserve p4est ka begin
        for i in 1:n_mirrors
            pq = pw_mirror_quadrant(pp, gp, i)
            dp = PointerWrapper(P4estPsData, pq.p.user_data[])
            ps_data = unsafe_pointer_to_objref(pointer(dp.ps_data))
            flag = if isa(ps_data, InsideSolidData) || ps_data.bound_enc < 0
                0.0
            else
                maximum(ps_data.lohner) > threshold ? 1.0 : 0.0
            end
            mirror_bufs[i] = Float64[flag]
        end
    end

    mirror_szs = fill(1, n_mirrors)
    ghost_szs  = fill(1, n_ghosts)
    recv_buf, recv_offsets =
        amr_exchange_varsize(p4est, ghost, mirror_bufs, mirror_szs, ghost_szs, 54)

    ghost_wrap = ka.kdata.ghost.ghost_wrap
    for i in 1:n_ghosts
        gw = ghost_wrap[i]
        isa(gw, GhostInsideSolidData) && continue
        gw.vs_data.vs_num == 0 && continue
        ghost_flags[objectid(gw)] = recv_buf[recv_offsets[i] + 1] > 0.5
    end
    return ghost_flags
end

# ---------------------------------------------------------------------------
# Ghost layer rebuild
# ---------------------------------------------------------------------------

function update_ghost!(p4est::Ptr{p4est_t}, ka::KA)
    kinfo = ka.kinfo
    if kinfo.forest.mesh != Ptr{p4est_mesh_t}(C_NULL)
        p4est_mesh_destroy(kinfo.forest.mesh)
        kinfo.forest.mesh = Ptr{p4est_mesh_t}(C_NULL)
    end
    if kinfo.forest.ghost != Ptr{p4est_ghost_t}(C_NULL)
        p4est_ghost_destroy(kinfo.forest.ghost)
    end
    kinfo.forest.ghost = p4est_ghost_new(p4est, P4EST_CONNECT_FULL)
    gb, gi = initialize_ghost_pool(p4est, kinfo)
    ka.kdata.ghost.ghost_buffer = gb
    ka.kdata.ghost.ghost_info   = gi
    ka.kdata.ghost.ghost_wrap   = initialize_ghost_wrap(kinfo, gb, gi)
end
function update_ghost!(p4est::Ptr{p8est_t}, ka::KA)
    kinfo = ka.kinfo
    if kinfo.forest.mesh != Ptr{p8est_mesh_t}(C_NULL)
        p8est_mesh_destroy(kinfo.forest.mesh)
        kinfo.forest.mesh = Ptr{p8est_mesh_t}(C_NULL)
    end
    if kinfo.forest.ghost != Ptr{p8est_ghost_t}(C_NULL)
        p8est_ghost_destroy(kinfo.forest.ghost)
    end
    kinfo.forest.ghost = p8est_ghost_new(p4est, P8EST_CONNECT_FULL)
    gb, gi = initialize_ghost_pool(p4est, kinfo)
    ka.kdata.ghost.ghost_buffer = gb
    ka.kdata.ghost.ghost_info   = gi
    ka.kdata.ghost.ghost_wrap   = initialize_ghost_wrap(kinfo, gb, gi)
end
