#!/bin/bash
set -euo pipefail

# -----------------------------------------------------------------------------
# hlpr_run_bamtofastq.sh
#
# Run parallel Cell Ranger bamtofastq using a Python helper script.
# Prefers execution inside the UV container, with fallback to host Python
# if Singularity is unavailable inside the container.
# Uses multiprocessing under the hood, --jobs argument allows to tweak number of workers.
# If none are provided, uses all available.
#
# Usage:
#   ./hlpr_run_bamtofastq.sh <DATASET_DIR> [--jobs N]
#
# Example:
#   ./hlpr_run_bamtofastq.sh SC2-23 --jobs 8
# -----------------------------------------------------------------------------

# ----------------------------- Args -------------------------------------------

if [ $# -lt 1 ]; then
    echo "Usage: $0 <DATASET_DIR> [--jobs N]"
    exit 1
fi

DATASET_DIR="$1"
shift
EXTRA_ARGS="$@"

# ----------------------------- Config -----------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

set -a
source "$SCRIPT_DIR/../../config.env"
set +a

# ----------------------------- Paths ------------------------------------------

PY_SCRIPT="../../resources/run_bamtofastq.py"
BAM_DIR="${DATADIR}/1a_source/${DATASET_DIR}"

if [ ! -f "$PY_SCRIPT" ]; then
    echo "ERROR: Python script not found: $PY_SCRIPT"
    exit 1
fi

if [ ! -d "$BAM_DIR" ]; then
    echo "ERROR: BAM directory not found: $BAM_DIR"
    exit 1
fi

# ----------------------------- Runtime Detection ------------------------------

echo "Checking Singularity availability inside UV container..."

if singularity exec "$UV_PATH" which singularity >/dev/null 2>&1; then
    RUN_MODE="uv"
    echo "Singularity found inside UV container"
else
    RUN_MODE="host"
    echo "Singularity NOT found inside UV container — falling back to host Python"
fi

# ----------------------------- Host Python Fallback ----------------------------

if [ "$RUN_MODE" = "host" ]; then
    if command -v python >/dev/null 2>&1; then
        HOST_PYTHON=python
    elif command -v python3 >/dev/null 2>&1; then
        HOST_PYTHON=python3
    else
        echo "ERROR: No python or python3 found on host"
        exit 1
    fi
    echo "Using host Python: $HOST_PYTHON"
fi

# ----------------------------- Execution --------------------------------------

echo "Running bamtofastq for dataset: $DATASET_DIR"
echo "BAM_DIR: $BAM_DIR"
echo "CELLRANGER_PATH: $CELLRANGER_PATH"
echo "Execution mode: $RUN_MODE"

if [ "$RUN_MODE" = "uv" ]; then
    singularity exec "$UV_PATH" \
        uv run python "$PY_SCRIPT" \
        "$CELLRANGER_PATH" \
        "$BAM_DIR" \
        $EXTRA_ARGS
else
    "$HOST_PYTHON" "$PY_SCRIPT" \
        "$CELLRANGER_PATH" \
        "$BAM_DIR" \
        $EXTRA_ARGS
fi

echo "bamtofastq processing completed."
