#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd)"
cd "$ROOT"

CORES="${RP3D_PROBE_CORES:-2048}"
TASKS_PER_NODE="${RP3D_TASKS_PER_NODE:-64}"
NODES="${RP3D_PROBE_NODES:-$(( (CORES + TASKS_PER_NODE - 1) / TASKS_PER_NODE ))}"
SBATCH_EXTRA="${RP3D_SBATCH_EXTRA:-}"

CASE_DIR="example/Riemann_problem/rp3d_scaling"
LOG_DIR="$CASE_DIR/logs"
RUN_SCRIPT="$CASE_DIR/slurm/run_case.sbatch"
mkdir -p "$LOG_DIR"

echo "submit probe cores=$CORES nodes=$NODES tasks_per_node=$TASKS_PER_NODE logs=$LOG_DIR"

# shellcheck disable=SC2086
sbatch $SBATCH_EXTRA \
    --job-name="rp3d_probe_${CORES}" \
    --output="$LOG_DIR/%x_%j.out" \
    --error="$LOG_DIR/%x_%j.err" \
    --nodes="$NODES" \
    --ntasks="$CORES" \
    --ntasks-per-node="$TASKS_PER_NODE" \
    "$RUN_SCRIPT" probe "$CORES"
