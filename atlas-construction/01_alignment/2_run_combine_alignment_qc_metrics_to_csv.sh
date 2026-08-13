#!/bin/bash
set -euo pipefail

# -----------------------------------------------------------------------------
# 2_run_cellranger_qc_metrics_to_csv.sh
# -----------------------------------------------------------------------------
# Runs cellranger_qc_metrics_to_csv.py inside a Singularity container using uv.
#
# Required config variables (from ../../config.env):
#   DATADIR
#   UV_PATH
# -----------------------------------------------------------------------------

############################
# Resolve paths
############################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY_SCRIPT="${SCRIPT_DIR}/../../resources/combine_alignment_qc_metrics_to_csv.py"

############################
# Load configuration
############################

set -a
source "$SCRIPT_DIR/../../config.env"
set +a

for var in DATADIR UV_PATH; do
    if [ -z "${!var:-}" ]; then
        echo "$(date): ERROR: $var not set in config.env"
        exit 1
    fi
done

############################
# Validate inputs
############################

INPUT_DIR="${DATADIR}/1c_aligned"

if [ ! -d "$INPUT_DIR" ]; then
    echo "$(date): ERROR: Input directory not found: $INPUT_DIR"
    exit 2
fi

if [ ! -f "$PY_SCRIPT" ]; then
    echo "$(date): ERROR: Python script not found: $PY_SCRIPT"
    exit 3
fi

if [ ! -f "$UV_PATH" ]; then
    echo "$(date): ERROR: Singularity image not found: $UV_PATH"
    exit 4
fi

############################
# Run QC aggregation
############################

echo "$(date): Running Cellranger QC metrics aggregation"
echo "$(date): Input directory: $INPUT_DIR"
echo "$(date): Python script: $PY_SCRIPT"
echo "$(date): Singularity image: $UV_PATH"

singularity exec "$UV_PATH" \
    uv run "$PY_SCRIPT" "$INPUT_DIR"

echo "$(date): QC metrics aggregation completed successfully"
