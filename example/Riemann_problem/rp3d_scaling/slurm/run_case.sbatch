#!/usr/bin/env bash
#SBATCH --job-name=rp3d
#SBATCH --exclusive
#SBATCH --time=02:00:00
#SBATCH --output=example/Riemann_problem/rp3d_scaling/logs/%x_%j.out
#SBATCH --error=example/Riemann_problem/rp3d_scaling/logs/%x_%j.err
##SBATCH --partition=your_partition
##SBATCH --account=your_account
#SBATCH --mem=0

set -euo pipefail

SCALING="${1:?usage: sbatch [slurm options] run_case.sbatch strong|weak|probe CORES}"
CORES="${2:?usage: sbatch [slurm options] run_case.sbatch strong|weak|probe CORES}"

if [[ "$SCALING" != "strong" && "$SCALING" != "weak" && "$SCALING" != "probe" ]]; then
    echo "SCALING must be 'strong', 'weak', or 'probe', got '$SCALING'" >&2
    exit 2
fi

ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "$ROOT"

# Optionally load modules without editing the script:
#   RP3D_MODULES="julia openmpi" sbatch ...
if [[ -n "${RP3D_MODULES:-}" ]]; then
    module purge
    for module_name in ${RP3D_MODULES}; do
        module load "$module_name"
    done
fi

if [[ "$SCALING" == "probe" ]]; then
    SCRIPT="example/Riemann_problem/rp3d_scaling/quick_probe_${CORES}.jl"
else
    SCRIPT="example/Riemann_problem/rp3d_scaling/${SCALING}_${CORES}.jl"
fi
SYSIMAGE="${KITAMR_SYSIMAGE:-$ROOT/example/Riemann_problem/rp3d_scaling/kitamr_rp3d_sys.so}"

if [[ ! -f "$SCRIPT" ]]; then
    echo "case script not found: $SCRIPT" >&2
    exit 2
fi

export JULIA_NUM_THREADS="${JULIA_NUM_THREADS:-1}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export JULIA_PKG_PRECOMPILE_AUTO="${JULIA_PKG_PRECOMPILE_AUTO:-0}"
export JULIA_PROJECT="${JULIA_PROJECT:-$ROOT}"

JULIA_ARGS=(--project="$JULIA_PROJECT")
if [[ -f "$SYSIMAGE" ]]; then
    JULIA_ARGS+=(-J "$SYSIMAGE")
else
    if [[ "${RP3D_ALLOW_NO_SYSIMAGE:-0}" == "1" ]]; then
        echo "warning: sysimage not found, running without -J: $SYSIMAGE" >&2
    else
        echo "sysimage not found: $SYSIMAGE" >&2
        echo "build it manually first, or set RP3D_ALLOW_NO_SYSIMAGE=1 for a smoke test without -J" >&2
        exit 2
    fi
fi

allocated_tasks="${SLURM_NTASKS:-$CORES}"
if [[ "$allocated_tasks" != "$CORES" ]]; then
    echo "warning: case expects $CORES MPI ranks, allocation has $allocated_tasks tasks" >&2
fi

echo "[$(date)] starting RP3D case"
echo "ROOT=$ROOT"
echo "SCALING=$SCALING CORES=$CORES"
echo "SCRIPT=$SCRIPT"
echo "SYSIMAGE=$SYSIMAGE"
echo "SLURM_JOB_ID=${SLURM_JOB_ID:-none}"
echo "SLURM_JOB_NUM_NODES=${SLURM_JOB_NUM_NODES:-unknown}"
echo "SLURM_NTASKS=${SLURM_NTASKS:-unknown}"
echo "SLURM_NTASKS_PER_NODE=${SLURM_NTASKS_PER_NODE:-unknown}"
echo "SLURM_JOB_NODELIST=${SLURM_JOB_NODELIST:-unknown}"
echo "SLURM_MEM_PER_NODE=${SLURM_MEM_PER_NODE:-unknown}"
echo "JULIA_NUM_THREADS=$JULIA_NUM_THREADS"
echo "RP3D_SETUP_TIME=${RP3D_SETUP_TIME:-default}"
echo "RP3D_SETUP_STEPS=${RP3D_SETUP_STEPS:-default}"
echo "RP3D_WARMUP_STEPS=${RP3D_WARMUP_STEPS:-default}"
echo "RP3D_TIMED_STEPS=${RP3D_TIMED_STEPS:-default}"
echo "RP3D_GC_BEFORE_TIMING=${RP3D_GC_BEFORE_TIMING:-default}"
echo "RP3D_GC_PASSES=${RP3D_GC_PASSES:-default}"
echo "RP3D_QUADRATURE_BOUND=${RP3D_QUADRATURE_BOUND:-default}"
echo "RP3D_TMAX_ESTIMATE=${RP3D_TMAX_ESTIMATE:-default}"

launcher="${RP3D_LAUNCHER:-srun}"
case "$launcher" in
    srun)
        SRUN_ARGS=()
        if [[ -n "${RP3D_SRUN_MPI:-}" ]]; then
            SRUN_ARGS+=(--mpi="$RP3D_SRUN_MPI")
        fi
        srun "${SRUN_ARGS[@]}" "${JULIA:-julia}" "${JULIA_ARGS[@]}" "$SCRIPT"
        ;;
    mpiexec)
        mpiexec -n "$CORES" "${JULIA:-julia}" "${JULIA_ARGS[@]}" "$SCRIPT"
        ;;
    mpirun)
        mpirun -n "$CORES" "${JULIA:-julia}" "${JULIA_ARGS[@]}" "$SCRIPT"
        ;;
    *)
        echo "unknown RP3D_LAUNCHER='$launcher' (use srun, mpiexec, or mpirun)" >&2
        exit 2
        ;;
esac

echo "[$(date)] RP3D case done"
