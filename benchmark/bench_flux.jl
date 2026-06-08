# =============================================================================
# bench_flux.jl  —  Quantitative baseline for the type-instability cost in flux!
# -----------------------------------------------------------------------------
# WHAT THIS MEASURES
#   The hot kernel `flux!(ka)` loops over every face of the local mesh and calls
#   `flux!(face, ka)` -> make_face_vs -> calc_flux/calc_domain_flux -> update_flux!.
#   The type-stability diagnosis (diag_typestab.jl) showed this chain is dynamic:
#     - `solver.flux` infers to `Type{<:AbstractFluxType}` (abstract)  -> RC-1
#     - `ka.kinfo` infers to `KInfo` (UnionAll, params dropped)        -> RC-2
#   => calc_flux returns `Tuple{Any,Any}`, so every face does a runtime dispatch
#      and BOXES its results. Boxing shows up as heap ALLOCATIONS that scale with
#      the number of faces. Therefore the cleanest, hardware-independent signal of
#      the instability is: *allocation count and bytes per `flux!(ka)` call*.
#
# WHAT WE REPORT (per MPI rank, then aggregated)
#   - #faces (and breakdown by face type) processed by `flux!(ka)`
#   - median & minimum wall time of `flux!(ka)`
#   - allocation COUNT and BYTES per call         <-- primary instability metric
#   - derived: allocs/face and bytes/face         <-- normalized, mesh-independent
#
# WHY flux!(ka) and not flux!(p4est,ka):
#   The (p4est,ka) variant overlaps MPI ghost/solid exchange; we want the pure
#   local compute kernel so the number isolates type instability, not comms.
#
# MESH IS FIXED & DETERMINISTIC so the baseline run and the post-fix run are
# directly comparable (identical face set). Run the SAME script before/after fix.
#
# USAGE
#   mpirun -np 1  julia --project=. bench_flux.jl baseline
#   mpirun -np 64 julia --project=. bench_flux.jl baseline      # under load
#   (after fixing RC-1+RC-2) re-run with tag `fixed`
# =============================================================================
using KitAMR, MPI
using BenchmarkTools
using Printf
MPI.Init()

const RANK  = MPI.Comm_rank(MPI.COMM_WORLD)
const NRANK = MPI.Comm_size(MPI.COMM_WORLD)
const TAG   = isempty(ARGS) ? "baseline" : ARGS[1]

# ---- deterministic mesh: 24x24 base trees, refine the left half twice ----
# (mixed FullFace / HangingFace / BackHangingFace / DomainFace, like real runs)
refine_left(midpoint, ds, kinfo, level) = (midpoint[1] < 0.0 && level < 2)

solver = Solver(;
    DIM = 2, NDF = 2,
    AMR_PS_MAXLEVEL = 2,
    AMR_VS_MAXLEVEL = 0,
    PS_DYNAMIC_AMR = false,
    VS_DYNAMIC_AMR = false,
    flux = CAIDVM,
    time_marching = CAIDVM_Marching,
)
gas = Gas(; K = 0.0, Kn = 0.075, ω = 0.81, ωᵣ = 0.81)
output = Output(solver)
udf = UDF(; static_ps_refine_flag = refine_left)
config = Configure(solver;
    geometry = [-0.5,0.5,-0.5,0.5],
    trees_num = [24,24],
    quadrature = [-5.0,5.0,-5.0,5.0],
    vs_trees_num = [16,16],
    IC = Uniform([1.,0.,0.,1.]),
    domain = [
        Domain(Maxwellian,1,[1.,0.,0.,1.]),
        Domain(Maxwellian,2,[1.,1.0*sqrt(5/6),0.,1.0]),
        Domain(Period,3), Domain(Period,4),
    ],
    output = output, gas = gas, user_defined = udf,
)

p4est, ka = initialize(config)

# ---- warm-up: populate slopes + JIT-compile the flux path (exclude compile time)
KitAMR.slope!(p4est, ka)
KitAMR.flux!(ka)              # force compilation of all face-type methods
KitAMR.slope!(p4est, ka)

# ---- face inventory (local to this rank) ----
faces = ka.kdata.field.faces
nfaces = length(faces)
tally = Dict{Symbol,Int}()
for f in faces
    k = nameof(typeof(f))
    tally[k] = get(tally, k, 0) + 1
end

# ---- benchmark the pure local kernel ----
# evals=1: one flux!(ka) per measured sample (we want per-call alloc/time).
# A generous sample budget gives a stable minimum/median.
b = @benchmark KitAMR.flux!($ka) samples=400 evals=1 seconds=20

t_med  = median(b).time      # ns
t_min  = minimum(b).time     # ns
allocs = b.allocs            # allocation COUNT per call
membytes = b.memory          # allocated BYTES per call

# ---- per-rank report ----
println("="^78)
@printf("[rank %d/%d] TAG=%s  nfaces=%d  %s\n", RANK, NRANK, TAG, nfaces, string(tally))
@printf("[rank %d] flux!(ka): t_min=%.1f µs  t_med=%.1f µs  allocs=%d  mem=%.1f KiB\n",
        RANK, t_min/1e3, t_med/1e3, allocs, membytes/1024)
@printf("[rank %d] normalized: %.2f allocs/face   %.1f bytes/face\n",
        RANK, allocs/max(nfaces,1), membytes/max(nfaces,1))

# ---- aggregate across ranks ----
tot_faces  = MPI.Allreduce(nfaces,   +, MPI.COMM_WORLD)
tot_allocs = MPI.Allreduce(allocs,   +, MPI.COMM_WORLD)
tot_mem    = MPI.Allreduce(membytes, +, MPI.COMM_WORLD)
max_tmed   = MPI.Allreduce(t_med,    max, MPI.COMM_WORLD)

if RANK == 0
    println("#"^78)
    @printf("# AGGREGATE TAG=%s over %d ranks:\n", TAG, NRANK)
    @printf("#   total faces/step      = %d\n", tot_faces)
    @printf("#   total allocs/step     = %d        (sum over ranks)\n", tot_allocs)
    @printf("#   total mem/step        = %.2f MiB  (sum over ranks)\n", tot_mem/1024/1024)
    @printf("#   slowest rank t_med    = %.1f µs\n", max_tmed/1e3)
    @printf("#   normalized            = %.2f allocs/face,  %.1f bytes/face\n",
            tot_allocs/max(tot_faces,1), tot_mem/max(tot_faces,1))
    println("#"^78)
    # append one machine-readable line for before/after comparison
    open("bench_flux_results.txt", "a") do io
        @printf(io, "%s\tnranks=%d\ttot_faces=%d\ttot_allocs=%d\ttot_mem_MiB=%.3f\tmax_tmed_us=%.2f\tallocs_per_face=%.3f\tbytes_per_face=%.2f\n",
                TAG, NRANK, tot_faces, tot_allocs, tot_mem, max_tmed/1e3,
                tot_allocs/max(tot_faces,1), tot_mem/max(tot_faces,1))
    end
end

KitAMR.finalize!(p4est, ka)
MPI.Finalize()
