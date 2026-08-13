#!/bin/bash
set -euo pipefail

###############################################################################
# Script Name: hlpr_gzip_files.sh
#
# Description:
#   Wrapper to run gzip_files.py inside the UV container.
#   Compresses all .fastq files in the directory:
#       ${DATADIR}/1a_source/${DATASET_DIR}
#
# Usage:
#   ./hlpr_gzip_files.sh <DATASET_DIR>
#
# Notes:
#   - Reads DATADIR and UV_PATH from ../../config.env
#   - Compresses all .fastq files that are not already gzipped
#   - Logs actions to <directory>/logs/compression.log
###############################################################################

# ----------------------------- Args -------------------------------------------

if [ $# -lt 1 ]; then
  echo "Usage: $0 <DATASET_DIR>"
  exit 1
fi

DATASET_DIR="$1"

# ----------------------------- Config -----------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

set -a
source "$SCRIPT_DIR/../../config.env"
set +a

if [ -z "$DATADIR" ]; then
  echo "ERROR: DATADIR not set in config.env"
  exit 1
fi

if [ -z "$UV_PATH" ]; then
  echo "ERROR: UV_PATH not set in config.env"
  exit 1
fi

# ----------------------------- Paths ------------------------------------------

TARGET_DIR="${DATADIR}/1a_source/${DATASET_DIR}"

if [ ! -d "$TARGET_DIR" ]; then
  echo "ERROR: Target directory '$TARGET_DIR' does not exist."
  exit 1
fi

SCRIPT_PATH="../../resources/gzip_files.py"

# ----------------------------- Execution --------------------------------------

echo "Running gzip_files.py on $TARGET_DIR using UV container at $UV_PATH"

singularity exec "$UV_PATH" python3 "$SCRIPT_PATH" "$TARGET_DIR"

echo "Done compressing .fastq files."
