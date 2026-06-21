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

mkdir -p "$LOG_DIR" "$OUTDIR"

mode_np() {
    case "$1" in
        both_on)  echo "${NP_BOTH_ON:-64}" ;;
        ps_only)  echo "${NP_PS_ONLY:-128}" ;;
        vs_only)  echo "${NP_VS_ONLY:-64}" ;;
        both_off) echo "${NP_BOTH_OFF:-128}" ;;
        *) echo "unknown mode $1" >&2; exit 2 ;;
    esac
}

for mode in both_on ps_only vs_only both_off; do
    np="$(mode_np "$mode")"
    nodes=$(( (np + TASKS_PER_NODE - 1) / TASKS_PER_NODE ))
    script="$CASE_DIR/${mode}.jl"
    echo "submit mode=$mode np=$np nodes=$nodes tasks_per_node=$TASKS_PER_NODE outdir=$OUTDIR"
    # shellcheck disable=SC2086
    sbatch $SBATCH_EXTRA \
        --job-name="r2d_vj_${mode}" \
        --output="$LOG_DIR/%x_%j.out" \
        --error="$LOG_DIR/%x_%j.err" \
        --nodes="$nodes" \
        --ntasks="$np" \
        --ntasks-per-node="$TASKS_PER_NODE" \
        --exclusive \
        "$RUNNER" "$script" VJE_OUTDIR="$OUTDIR" VJE_LABEL="$mode"
done
