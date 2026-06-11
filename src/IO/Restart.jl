# Restart / checkpoint.
#
# `save_for_restart` writes the *full phase-space* state (every velocity cell's distribution `df`,
# the forest topology, the configuration and the dynamic status) so that a run can be resumed
# later by `restart` — possibly on a *different* number of MPI ranks. This is heavier than
# [`save_result`](@ref) (which only stores macroscopic fields + user-selected velocity spaces for
# post-processing) and is the only output from which a simulation can continue bit-for-bit.
#
# On-disk layout (under `dir_path/`, default `./restart/`):
#   p.*                forest topology + global quadrant ordering (p4est_save_ext, save_data = 0)
#   solverset.jld2     SolverSet(ConfigureForSave(config), mpi_size)  — config incl. IC/UDF
#   output.jld2        the Output struct (not carried by ConfigureForSave)
#   status.jld2        dynamic Status scalars needed to resume the time loop
#   manifest.jld2      format/version + the gfq partition map + per-rank counts (integrity)
#   restart_<rank>.jld2  per-rank, compressed: per-cell phase space (vs_nums/bound_encs/ws + vs df)

const RESTART_FORMAT_VERSION = 1

# ---------------------------------------------------------------------------------------------------
# Small MPI helpers
# ---------------------------------------------------------------------------------------------------
function _bcast_bool(b::Bool)
    a = Cint[b]
    MPI.Bcast!(a, 0, MPI.COMM_WORLD)
    return a[1] != 0
end
# Broadcast a String from rank 0 to every rank (used to raise the same integrity error everywhere).
function _bcast_string(s::String)
    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        v = collect(s)
        MPI.Bcast!(Int[length(v)], 0, MPI.COMM_WORLD)
        MPI.Bcast!(v, 0, MPI.COMM_WORLD)
        return s
    else
        n = Int[0]
        MPI.Bcast!(n, 0, MPI.COMM_WORLD)
        v = Vector{Char}(undef, n[1])
        MPI.Bcast!(v, 0, MPI.COMM_WORLD)
        return String(v)
    end
end
function _human_bytes(n::Real)
    units = ("B", "KiB", "MiB", "GiB", "TiB")
    x = float(n); i = 1
    while x >= 1024 && i < length(units)
        x /= 1024; i += 1
    end
    return string(round(x; digits = 2), " ", units[i])
end

# ---------------------------------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------------------------------
"""
$(TYPEDSIGNATURES)
Write a full phase-space checkpoint of the current state to `dir_path` (default `./restart`), from
which the run can be resumed with [`restart`](@ref) — including on a different number of MPI ranks.

Unlike [`save_result`](@ref), this stores every velocity cell's distribution function, so the
checkpoint can be large. Each rank estimates its byte footprint; rank 0 prints the aggregate size
and, when `confirm = true` **and** `stdin` is an interactive terminal, asks for `y/n` confirmation
before writing (non-interactive runs just print the size and proceed). Returns the directory path,
or `nothing` if the user declines.

Files written (see the module comment for the full layout): the p4est forest (`p.*`), the
configuration (`solverset.jld2` incl. the initial-condition / user-defined functions, plus
`output.jld2`), the dynamic status (`status.jld2`), an integrity `manifest.jld2`, and one compressed
`restart_<rank>.jld2` per rank carrying the per-cell phase space.

Call after a time step (e.g. inside or after [`solve!`](@ref)); the saved distribution `df` is taken
as-is, so no slope/macro update is performed.
"""
function save_for_restart(p4est::P_pxest_t, ka::KA{DIM,NDF}; dir_path::String = "restart", confirm::Bool = true) where{DIM,NDF}
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    mpi_size = MPI.Comm_size(comm)
    fp = PointerWrapper(p4est)
    N = fp.local_num_quadrants[]
    trees = ka.kdata.field.trees.data

    # First pass: per-cell scalars + total velocity-cell count.
    vs_nums = Vector{Int}(undef, N)
    bound_encs = Vector{Int}(undef, N)
    ws = zeros(Float64, N, DIM + 2)
    index = 1
    for tree in trees
        for ps_data in tree
            if isa(ps_data, InsideSolidData)
                vs_nums[index] = 0
                bound_encs[index] = ps_data.bound_enc
            else
                vs_nums[index] = ps_data.vs_data.vs_num
                bound_encs[index] = ps_data.bound_enc
                ws[index, :] .= ps_data.w
            end
            index += 1
        end
    end
    Nv = sum(vs_nums)

    # Second pass: flattened velocity-space structure + distribution.
    vs_levels = Vector{Int8}(undef, Nv)
    vs_midpoints = Matrix{Float64}(undef, Nv, DIM)
    vs_df = Matrix{Float64}(undef, Nv, NDF)
    offset = 0
    for tree in trees
        for ps_data in tree
            isa(ps_data, InsideSolidData) && continue
            vs_data = ps_data.vs_data
            vn = vs_data.vs_num
            rng = offset+1:offset+vn
            vs_levels[rng] .= vs_data.level
            vs_midpoints[rng, :] .= vs_data.midpoint
            vs_df[rng, :] .= vs_data.df
            offset += vn
        end
    end

    # Size estimate (raw, uncompressed) and optional confirmation.
    local_bytes = sizeof(vs_nums) + sizeof(bound_encs) + sizeof(ws) +
                  sizeof(vs_levels) + sizeof(vs_midpoints) + sizeof(vs_df)
    total_bytes = MPI.Reduce(local_bytes, +, 0, comm)
    proceed = true
    if rank == 0
        println("Restart checkpoint: ≈ $(_human_bytes(total_bytes)) of phase-space data " *
                "($(mpi_size) rank file(s)) → $(dir_path)/")
        if confirm && isa(stdin, Base.TTY)
            print("Continue saving? [y/N]: ")
            ans = strip(readline())
            proceed = (ans == "y" || ans == "Y" || ans == "yes")
        end
    end
    proceed = _bcast_bool(proceed)
    if !proceed
        rank == 0 && println("Restart save aborted by user.")
        return nothing
    end

    # Gather per-rank counts for the manifest.
    local_quad_nums = MPI.Gather(N, 0, comm)
    vs_totals = MPI.Gather(Nv, 0, comm)

    dir = endswith(dir_path, "/") ? dir_path : dir_path * "/"
    if rank == 0
        !isdir(dir) && mkpath(dir)
    end
    MPI.Barrier(comm)

    # Save the forest (topology only; per-quad user_data is a meaningless Julia pointer).
    # Absolute path instead of cd("p"): cd mutates the process-global CWD (unsafe and
    # not restored if the collective save throws); the Cstring filename takes a full path.
    fpath = joinpath(abspath(dir), "p")
    if DIM == 2
        GC.@preserve p4est p4est_save_ext(fpath, p4est, Cint(0), Cint(0))
    else
        GC.@preserve p4est p8est_save_ext(fpath, p4est, Cint(0), Cint(0))
    end

    # Rank-0 metadata.
    if rank == 0
        config = ka.kinfo.config
        status = ka.kinfo.status
        gfq = unsafe_wrap(Vector{Int}, pointer(fp.global_first_quadrant), mpi_size + 1)
        jldopen(dir * "solverset.jld2", "w") do f
            f["solverset"] = SolverSet(ConfigureForSave(config), mpi_size)
        end
        jldopen(dir * "output.jld2", "w") do f
            f["output"] = config.output
        end
        jldopen(dir * "status.jld2", "w") do f
            f["step"] = status.step
            f["sim_time"] = status.sim_time
            f["Dt"] = status.Δt
            f["Dt_xi"] = status.Δt_ξ
            f["gradmax"] = status.gradmax
            f["ps_adapt_step"] = status.ps_adapt_step
            f["vs_adapt_step"] = status.vs_adapt_step
            f["partition_step"] = status.partition_step
        end
        jldopen(dir * "manifest.jld2", "w") do f
            f["format_version"] = RESTART_FORMAT_VERSION
            f["DIM"] = DIM
            f["NDF"] = NDF
            f["mpi_size"] = mpi_size
            f["global_quad_num"] = fp.global_num_quadrants[]
            f["gfq"] = collect(gfq)
            f["local_quad_nums"] = local_quad_nums
            f["vs_totals"] = vs_totals
        end
    end

    # Per-rank phase space (compressed).
    jldopen(dir * "restart_" * string(rank) * ".jld2", "w"; compress = true) do f
        f["vs_nums"] = vs_nums
        f["bound_encs"] = bound_encs
        f["ws"] = ws
        f["vs_levels"] = vs_levels
        f["vs_midpoints"] = vs_midpoints
        f["vs_df"] = vs_df
    end
    MPI.Barrier(comm)
    rank == 0 && println("Restart checkpoint written to $(dir).")
    return dir
end

# ---------------------------------------------------------------------------------------------------
# Restart
# ---------------------------------------------------------------------------------------------------
# Load the forest saved by `save_for_restart` (topology only). `p4est_load_ext` allocates and fills the
# connectivity through `cnn`; it lives inside the returned forest and is freed by `finalize!`.
function _restart_load_forest(dir::String, DIM::Integer)
    # The forest was saved with save_data = 0 (no per-quad section). `*_load_ext` therefore restores a
    # forest whose quadrant data_size is 0 — so we load topology only (data_size = 0, load_data = 0),
    # then `*_reset_data` to allocate a fresh per-quad `user_data` layer of `sizeof(P4estPsData)` into
    # which `_restart_attach!` writes the Julia `P4estPsData` pointers. autopartition = 1 lets p4est
    # balance the load across the (possibly different) current number of ranks.
    # Absolute path (matches save_for_restart) instead of cd("p"): no process-global CWD mutation.
    fpath = joinpath(abspath(dir), "p")
    dsz = Csize_t(sizeof(P4estPsData))
    if DIM == 2
        cnn = Ptr{Ptr{p4est_connectivity_t}}(Libc.malloc(sizeof(Ptr{Ptr{p4est_connectivity_t}})))
        p4est = GC.@preserve cnn p4est_load_ext(fpath, MPI.COMM_WORLD, Cint(0), Cint(0), Cint(1), Cint(0), C_NULL, cnn)
        GC.@preserve p4est p4est_reset_data(p4est, dsz, C_NULL, C_NULL)
    else
        cnn = Ptr{Ptr{p8est_connectivity_t}}(Libc.malloc(sizeof(Ptr{Ptr{p8est_connectivity_t}})))
        p4est = GC.@preserve cnn p8est_load_ext(fpath, MPI.COMM_WORLD, Cint(0), Cint(0), Cint(1), Cint(0), C_NULL, cnn)
        GC.@preserve p4est p8est_reset_data(p4est, dsz, C_NULL, C_NULL)
    end
    return p4est
end

# Rank-0 integrity check; returns an error message ("" when everything is consistent).
function _restart_integrity_message(dir::String, manifest)
    problems = String[]
    old_size = manifest["mpi_size"]
    manifest["format_version"] == RESTART_FORMAT_VERSION ||
        push!(problems, "manifest format_version $(manifest["format_version"]) ≠ expected $(RESTART_FORMAT_VERSION)")
    for f in ("solverset.jld2", "output.jld2", "status.jld2")
        isfile(dir * f) || push!(problems, "missing metadata file: $(dir*f)")
    end
    isempty(filter(f -> startswith(f, "p"), readdir(dir))) &&
        push!(problems, "missing p4est forest files (prefix 'p') in $(dir)")
    lqn = manifest["local_quad_nums"]; vst = manifest["vs_totals"]
    for r in 0:old_size-1
        f = dir * "restart_" * string(r) * ".jld2"
        if !isfile(f)
            push!(problems, "missing per-rank data file: $f")
            continue
        end
        try
            vsn = jldopen(io -> io["vs_nums"], f)        # cheap per-cell array
            length(vsn) == lqn[r+1] ||
                push!(problems, "$f: cell count $(length(vsn)) ≠ manifest $(lqn[r+1])")
            sum(vsn) == vst[r+1] ||
                push!(problems, "$f: velocity-cell total $(sum(vsn)) ≠ manifest $(vst[r+1])")
        catch e
            push!(problems, "$f: failed to read ($e)")
        end
    end
    return isempty(problems) ? "" : "Restart integrity check failed:\n  - " * join(problems, "\n  - ")
end

# Reconstruct this rank's ordered ps_data list from the saved per-cell arrays, mapping the saved
# partition (gfq_old, old_size files) onto this rank's global quadrant range under the *current*
# communicator. p4est save/load preserves the global Morton ordering, so concatenating the
# overlapping old-rank files in order and slicing out [G0, G1) is exact.
function _restart_build_ps_list(dir::String, kinfo::KInfo{DIM,NDF}, gfq_old::Vector{Int}, p4est::P_pxest_t) where{DIM,NDF}
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    fp = PointerWrapper(p4est)
    gfq_new = unsafe_wrap(Vector{Int}, pointer(fp.global_first_quadrant), MPI.Comm_size(comm) + 1)
    G0 = gfq_new[rank+1]; G1 = gfq_new[rank+2]
    count = G1 - G0
    ps_list = Vector{AbstractPsData{DIM,NDF}}(undef, count)
    count == 0 && return ps_list

    first_old = searchsortedlast(gfq_old, G0) - 1          # old rank owning global cell G0
    last_old = searchsortedlast(gfq_old, G1 - 1) - 1       # old rank owning global cell G1-1

    vs_nums = Int[]; bound_encs = Int[]
    ws = Matrix{Float64}(undef, 0, DIM + 2)
    vs_levels = Int8[]
    vs_midpoints = Matrix{Float64}(undef, 0, DIM)
    vs_df = Matrix{Float64}(undef, 0, NDF)
    for r in first_old:last_old
        d = load(dir * "restart_" * string(r) * ".jld2")
        append!(vs_nums, d["vs_nums"]); append!(bound_encs, d["bound_encs"])
        ws = vcat(ws, d["ws"])
        append!(vs_levels, d["vs_levels"])
        vs_midpoints = vcat(vs_midpoints, d["vs_midpoints"])
        vs_df = vcat(vs_df, d["vs_df"])
    end

    base = gfq_old[first_old+1]          # global index of the first concatenated cell
    sa = G0 - base + 1                   # 1-based first cell for this rank
    # cumulative velocity-cell offsets over the concatenation
    cum = cumsum(vs_nums)
    vs_start = sa == 1 ? 0 : cum[sa-1]   # velocity cells before the first selected cell

    quadrature = kinfo.config.quadrature
    vs_space = 1.0
    for i = 1:DIM
        vs_space *= quadrature[2i] - quadrature[2i-1]
    end
    tree_weight = vs_space / reduce(*, kinfo.config.vs_trees_num)

    voff = vs_start
    for i in 1:count
        ci = sa + i - 1
        vn = vs_nums[ci]
        if vn == 0
            ps_list[i] = InsideSolidData(DIM, NDF; bound_enc = bound_encs[ci])
        else
            w = ws[ci, :]
            levels = vs_levels[voff+1:voff+vn]
            mids = vs_midpoints[voff+1:voff+vn, :]
            df = vs_df[voff+1:voff+vn, :]
            weight = @. tree_weight / 2.0^(DIM * levels)
            vs_data = VsData{DIM,NDF}(vn, levels, weight, mids, df, zeros(vn, NDF, DIM), zeros(vn, NDF))
            ps_list[i] = PsData(DIM, NDF; bound_enc = bound_encs[ci], w = w, prim = get_prim(w, kinfo), vs_data = vs_data)
        end
        voff += vn
    end
    return ps_list
end

# Attach the reconstructed ps_data to the loaded quadrants (in local order) and build the trees.
function _restart_attach!(p4est::P_pxest_t, kinfo::KInfo{DIM,NDF}, ps_list::Vector{AbstractPsData{DIM,NDF}}) where{DIM,NDF}
    fp = PointerWrapper(p4est)
    trees_data = [AbstractPsData{DIM,NDF}[] for _ in 1:fp.last_local_tree[]-fp.first_local_tree[]+1]
    trees = PsTrees{DIM,NDF}(trees_data, fp.first_local_tree[] - 1)
    ctx = Any[trees, ps_list, Ref(0)]
    p_data = pointer_from_objref(ctx)
    GC.@preserve ctx AMR_volume_iterate(p4est; user_data = p_data) do ip, data, dp
        trees, ps_list, cnt = unsafe_pointer_to_objref(data)
        cnt[] += 1
        ps_data = ps_list[cnt[]]
        ds, midpoint = quad_to_cell(ip.p4est, ip.treeid[], ip.quad)
        ps_data.ds = ds
        ps_data.midpoint = midpoint
        isa(ps_data, PsData) && (ps_data.quadid = global_quadid(ip))
        dp[] = P4estPsData(pointer_from_objref(ps_data))
        push!(trees.data[ip.treeid[]-trees.offset], ps_data)
    end
    return trees
end

"""
$(TYPEDSIGNATURES)
Resume a simulation from a checkpoint written by [`save_for_restart`](@ref) and return `(p4est, ka)`,
ready to hand to [`solve!`](@ref). The number of MPI ranks may differ from the run that wrote the
checkpoint; the saved field is remapped onto the current partition.

By default the configuration (including the initial-condition and user-defined functions) is loaded
from `dir_path/solverset.jld2`. **Those named functions must be defined in the current session**
(i.e. `include` the same user-defined-function script before calling `restart`), otherwise JLD2
cannot resolve them. Alternatively pass `config = <your live Configure>` to use the configuration
from your script verbatim; its `DIM`/`NDF`/`geometry`/`quadrature`/`vs_trees_num` are then checked
against the manifest.

`check_integrity = true` (default) verifies, on rank 0, that the manifest, configuration, forest and
all per-rank data files are present and mutually consistent, raising a descriptive error on every
rank if not.

The loaded distribution `df` is the restored state — no initial condition is re-applied and no
pre-refinement is performed. Reports the resumed status before returning.
"""
function restart(dir_path::String; config = nothing, check_integrity::Bool = true)
    !MPI.Initialized() && MPI.Init()
    comm = MPI.COMM_WORLD
    dir = endswith(dir_path, "/") ? dir_path : dir_path * "/"

    manifest = load(dir * "manifest.jld2")
    DIM = manifest["DIM"]; NDF = manifest["NDF"]

    if check_integrity
        msg = MPI.Comm_rank(comm) == 0 ? _restart_integrity_message(dir, manifest) : ""
        msg = _bcast_string(msg)
        isempty(msg) || error(msg)
    end

    # Configuration: from disk (self-contained) or user-supplied (overriding).
    if config === nothing
        solverset = load(dir * "solverset.jld2")["solverset"]
        cfs = solverset.config
        output = isfile(dir * "output.jld2") ? load(dir * "output.jld2")["output"] : Output(cfs.solver)
        config = Configure{DIM,NDF}(cfs.geometry, cfs.trees_num, cfs.quadrature, cfs.vs_trees_num,
            cfs.IC, cfs.domain, cfs.IB, cfs.gas, cfs.solver, output, cfs.user_defined)
    else
        (typeof(config).parameters[1] == DIM && typeof(config).parameters[2] == NDF) ||
            error("supplied `config` is $(typeof(config).parameters[1])D/NDF=$(typeof(config).parameters[2]); checkpoint is $(DIM)D/NDF=$NDF.")
    end

    kinfo = KInfo(config)
    p4est = _restart_load_forest(dir, DIM)
    kinfo.forest.p4est = p4est

    ps_list = _restart_build_ps_list(dir, kinfo, manifest["gfq"], p4est)
    trees = _restart_attach!(p4est, kinfo, ps_list)

    kdata = KData(trees)
    ka = KA(kinfo, kdata)
    PointerWrapper(p4est).user_pointer = pointer_from_objref(ka)

    # Same setup tail as `initialize`, minus the IC and the pre-refinement: the saved `df` is the state.
    initialize_forest!(p4est, kinfo)
    kdata.ghost = initialize_ghost(p4est, kinfo)
    initialize_neighbor_data!(p4est, ka)
    MPI.Comm_size(comm) > 1 && vs_ghost_exchange!(p4est, ka)
    update_neighbor!(p4est, ka)
    initialize_solid_neighbor!(ka)
    initialize_faces!(p4est, ka)
    initialize_immersed_boundaries!(ka)

    # Restore the dynamic status onto a config-consistent fresh Status.
    st = load(dir * "status.jld2")
    status = ka.kinfo.status
    status.step = st["step"]
    status.sim_time = st["sim_time"]
    status.Δt = st["Dt"]
    status.Δt_ξ = st["Dt_xi"]
    status.gradmax = st["gradmax"]
    status.ps_adapt_step = st["ps_adapt_step"]
    status.vs_adapt_step = st["vs_adapt_step"]
    status.partition_step = st["partition_step"]

    MPI.Comm_rank(comm) == 0 && println("Restarted from $(dir) at step $(status.step), sim_time $(status.sim_time).")
    execute_check!(p4est, ka)
    return p4est, ka
end
