#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
ROOT="$(cd "$CASE_DIR/../.." && pwd)"
cd "$ROOT"

RUNNER="${RUNNER:-$ROOT/slurm_templates/julia_mpi.sbatch}"
LOG_DIR="${LOG_DIR:-$CASE_DIR/logs}"
OUTDIR="${VJE_OUTDIR:-$CASE_DIR/out}"
SBATCH_EXTRA="${SBATCH_EXTRA:-}"

mkdir -p "$LOG_DIR"

echo "submit analysis outdir=$OUTDIR"

MPI_LAUNCHER=none \
    sbatch $SBATCH_EXTRA \
    --job-name="r2d_vj_analysis" \
    --output="$LOG_DIR/%x_%j.out" \
    --error="$LOG_DIR/%x_%j.err" \
    --nodes=1 \
    --ntasks=1 \
    "$RUNNER" "$CASE_DIR/analyze.jl" "$OUTDIR"
