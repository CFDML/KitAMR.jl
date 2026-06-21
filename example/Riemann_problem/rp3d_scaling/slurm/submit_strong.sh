#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd)"
cd "$ROOT"

CORES_LIST=("$@")
if [[ ${#CORES_LIST[@]} -eq 0 ]]; then
    CORES_LIST=(64 128 256 512 1024 2048)
fi

# Strong scaling keeps one large fixed problem.  Use enough nodes for memory,
# then place fewer MPI ranks per node for the smaller core counts.
NODES="${RP3D_STRONG_NODES:-32}"
CORES_PER_NODE="${RP3D_CORES_PER_NODE:-64}"
SBATCH_EXTRA="${RP3D_SBATCH_EXTRA:-}"
SCRIPT="example/Riemann_problem/rp3d_scaling/slurm/run_case.sbatch"
CASE_DIR="example/Riemann_problem/rp3d_scaling"
LOG_DIR="$CASE_DIR/logs"
mkdir -p "$LOG_DIR"

for cores in "${CORES_LIST[@]}"; do
    tasks_per_node=$(( (cores + NODES - 1) / NODES ))
    if (( tasks_per_node > CORES_PER_NODE )); then
        echo "cores=$cores needs tasks_per_node=$tasks_per_node > CORES_PER_NODE=$CORES_PER_NODE" >&2
        exit 2
    fi
    echo "submit strong cores=$cores nodes=$NODES tasks_per_node=$tasks_per_node"
    # shellcheck disable=SC2086
    sbatch $SBATCH_EXTRA \
        --job-name="rp3d_s${cores}" \
        --output="$LOG_DIR/%x_%j.out" \
        --error="$LOG_DIR/%x_%j.err" \
        --nodes="$NODES" \
        --ntasks="$cores" \
        --ntasks-per-node="$tasks_per_node" \
        "$SCRIPT" strong "$cores"
done
