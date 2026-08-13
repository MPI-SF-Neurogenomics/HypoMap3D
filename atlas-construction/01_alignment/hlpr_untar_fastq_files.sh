#!/bin/bash
set -euo pipefail

# -----------------------------------------------------------------------------
# hlpr_untar_fastq_files.sh
# -----------------------------------------------------------------------------
# Runs untar_fastq_files.py inside the UV Singularity container.
#
# Argument:
#   1) DATASET_DIR   (e.g. SC2-23)
#
# The script constructs:
#   INPUT_DIR = ${DATADIR}/1a_source/${DATASET_DIR}
#
# Required config variables (from ../../config.env):
#   DATADIR
#   UV_PATH
#
# Python script location:
#   ../../resources/untar_fastq_files.py
# -----------------------------------------------------------------------------

############################
# Arguments
############################

if [ $# -ne 1 ]; then
    echo "Usage: $0 <DATASET_DIR>"
    exit 1
fi

DATASET_DIR="$1"

############################
# Resolve paths
############################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY_SCRIPT="${SCRIPT_DIR}/../../resources/untar_fastq_files.py"

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
# Construct input directory
############################

INPUT_DIR="${DATADIR}/1a_source/${DATASET_DIR}"

if [ ! -d "$INPUT_DIR" ]; then
    echo "$(date): ERROR: Input directory not found: $INPUT_DIR"
    exit 1
fi

if [ ! -f "$PY_SCRIPT" ]; then
    echo "$(date): ERROR: Python script not found: $PY_SCRIPT"
    exit 1
fi

############################
# Run untar script
############################

echo "Running untar_fastq_files.py"
echo "  Dataset:   $DATASET_DIR"
echo "  Input dir: $INPUT_DIR"
echo "  Container: $UV_PATH"

singularity exec "$UV_PATH" \
    uv run python "$PY_SCRIPT" "$INPUT_DIR"

echo "untar_fastq_files completed successfully for dataset ${DATASET_DIR}"
