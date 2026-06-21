# =============================================================================
# build_sysimage.jl — Bake the RP3D scaling workload into a sysimage so per-rank
# JIT'd native code becomes a shared mmap'd .so instead of private memory
# replicated on every MPI rank.
#
# BUILD on a node with the same CPU architecture as the compute nodes:
#   julia -e 'using Pkg; Pkg.add("PackageCompiler")'
#   julia --project=. example/Riemann_problem/rp3d_scaling/build_sysimage.jl
#
# USE:
#   mpiexec -n 2048 julia --project=. -J example/Riemann_problem/rp3d_scaling/kitamr_rp3d_sys.so \
#       example/Riemann_problem/rp3d_scaling/strong_2048.jl
#
# NOTES: rebuild after updating KitAMR or its dependencies.  On heterogeneous
# systems, either build on the target compute-node microarchitecture or set
# RP3D_CPU_TARGET to a portable target accepted by Julia's --cpu-target.
# =============================================================================

using PackageCompiler

# The sysimage freezes the MPI binary choice.  On SLURM/srun clusters this should
# usually be "system"; otherwise the image can bake a JLL MPI and fail at startup.
import MPIPreferences
MPIPreferences.binary == "system" || @warn """
  MPIPreferences.binary = $(MPIPreferences.binary) (not "system").
  On an srun cluster, configure MPIPreferences.use_system_binary() in this
  project before building, or the sysimage may bake a JLL MPI that fails under
  the launcher used for production runs."""

const HERE = @__DIR__
const DEFAULT_SYSIMAGE = joinpath(HERE, "kitamr_rp3d_sys.so")
const SYSIMAGE = get(ENV, "KITAMR_SYSIMAGE", DEFAULT_SYSIMAGE)

# PackageCompiler executes precompile_execution_file in a single Julia process.
# That covers most kernel code, but not all Comm_size > 1 branches: p4est
# partition, ghost exchange, nonblocking MPI sends/receives, and global load
# reductions.  Generate a small two-rank trace once, then rebuild with it:
#
#   mpiexec -n 2 julia --project=. \
#     --trace-compile=example/Riemann_problem/rp3d_scaling/rp3d_stmts.jl \
#     example/Riemann_problem/rp3d_scaling/precompile_workload.jl
#
# On SLURM clusters, replace mpiexec by the production launcher, for example:
#
#   srun -n 2 julia --project=. \
#     --trace-compile=example/Riemann_problem/rp3d_scaling/rp3d_stmts.jl \
#     example/Riemann_problem/rp3d_scaling/precompile_workload.jl
#
# This script auto-detects rp3d_stmts.jl beside itself. Override with
# RP3D_PRECOMPILE_STMTS=/path/to/stmts.jl.
const STMTS = get(ENV, "RP3D_PRECOMPILE_STMTS", joinpath(HERE, "rp3d_stmts.jl"))
stmts_kw = isfile(STMTS) ? (; precompile_statements_file = STMTS) : (;)
isfile(STMTS) ? @info("including multi-rank statements from $STMTS") :
    @warn("no precompile_statements_file ($STMTS): MPI communication paths will JIT at run time; generate a 2-rank trace first for production runs")

build_args = String["--strip-metadata"]
if haskey(ENV, "RP3D_CPU_TARGET") && !isempty(ENV["RP3D_CPU_TARGET"])
    push!(build_args, "--cpu-target=$(ENV["RP3D_CPU_TARGET"])")
end

create_sysimage(
    [:KitAMR, :MPI];
    sysimage_path = SYSIMAGE,
    precompile_execution_file = joinpath(HERE, "precompile_workload.jl"),
    stmts_kw...,
    # Critical on system-MPI clusters: avoid baking the full dependency closure,
    # which may include OpenMPI_jll/MPICH_jll and initialize the wrong libmpi.
    include_transitive_dependencies = false,
    sysimage_build_args = Cmd(build_args),
)

println("\nWrote ", SYSIMAGE)
