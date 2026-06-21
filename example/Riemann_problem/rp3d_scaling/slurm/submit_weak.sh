#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd)"
cd "$ROOT"

CORES_LIST=("$@")
if [[ ${#CORES_LIST[@]} -eq 0 ]]; then
    CORES_LIST=(64 128 256 512 1024 2048)
fi

TASKS_PER_NODE="${RP3D_TASKS_PER_NODE:-64}"
SBATCH_EXTRA="${RP3D_SBATCH_EXTRA:-}"
SCRIPT="example/Riemann_problem/rp3d_scaling/slurm/run_case.sbatch"
CASE_DIR="example/Riemann_problem/rp3d_scaling"
LOG_DIR="$CASE_DIR/logs"
mkdir -p "$LOG_DIR"

for cores in "${CORES_LIST[@]}"; do
    nodes=$(( (cores + TASKS_PER_NODE - 1) / TASKS_PER_NODE ))
    echo "submit weak cores=$cores nodes=$nodes tasks_per_node=$TASKS_PER_NODE"
    # shellcheck disable=SC2086
    sbatch $SBATCH_EXTRA \
        --job-name="rp3d_w${cores}" \
        --output="$LOG_DIR/%x_%j.out" \
        --error="$LOG_DIR/%x_%j.err" \
        --nodes="$nodes" \
        --ntasks="$cores" \
        --ntasks-per-node="$TASKS_PER_NODE" \
        "$SCRIPT" weak "$cores"
done
