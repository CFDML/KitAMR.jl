#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
ROOT="$(cd "$CASE_DIR/../.." && pwd)"
cd "$ROOT"

RUNNER="${RUNNER:-$ROOT/slurm_templates/julia_mpi.sbatch}"
LOG_DIR="${LOG_DIR:-$CASE_DIR/logs}"
OUTDIR="${VJE_OUTDIR:-$CASE_DIR/out}"
TASKS_PER_NODE="${TASKS_PER_NODE:-64}"
SBATCH_EXTRA="${SBATCH_EXTRA:-}"
DIM="${VJE_DIM:-2}"

mkdir -p "$LOG_DIR" "$OUTDIR"

mode_np() {
    case "$1" in
        both_on)  echo "${NP_BOTH_ON:-64}" ;;
        both_off) echo "${NP_BOTH_OFF:-1024}" ;;
        both_on_3d)  echo "${NP_BOTH_ON_3D:-${NP_BOTH_ON:-64}}" ;;
        both_off_3d) echo "${NP_BOTH_OFF_3D:-${NP_BOTH_OFF:-1024}}" ;;
        *) echo "unknown mode $1" >&2; exit 2 ;;
    esac
}

if (( $# > 0 )); then
    MODES=("$@")
elif [[ -n "${VJE_MODES:-}" ]]; then
    read -r -a MODES <<< "${VJE_MODES//,/ }"
else
    case "$DIM" in
        2) MODES=(both_on both_off) ;;
        3) MODES=(both_on_3d) ;;
        *) echo "VJE_DIM must be 2 or 3, got $DIM" >&2; exit 2 ;;
    esac
fi

for mode in "${MODES[@]}"; do
    np="$(mode_np "$mode")"
    nodes=$(( (np + TASKS_PER_NODE - 1) / TASKS_PER_NODE ))
    script="$CASE_DIR/${mode}.jl"
    if [[ "$mode" == *_3d ]]; then
        job_name="r3d_vj_${mode%_3d}"
    else
        job_name="r2d_vj_${mode}"
    fi
    echo "submit mode=$mode dim=$DIM np=$np nodes=$nodes tasks_per_node=$TASKS_PER_NODE outdir=$OUTDIR"
    # shellcheck disable=SC2086
    sbatch $SBATCH_EXTRA \
        --constraint="64core" \
        --job-name="$job_name" \
        --output="$LOG_DIR/%x_%j.out" \
        --error="$LOG_DIR/%x_%j.err" \
        --nodes="$nodes" \
        --ntasks="$np" \
        --ntasks-per-node="$TASKS_PER_NODE" \
        --exclusive \
        --mem=0 \
        "$RUNNER" "$script" VJE_DIM="$DIM" VJE_OUTDIR="$OUTDIR" VJE_LABEL="$mode"
done
