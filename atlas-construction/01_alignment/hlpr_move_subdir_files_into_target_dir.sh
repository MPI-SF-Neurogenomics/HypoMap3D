#!/bin/bash
set -euo pipefail

###############################################################################
# Script Name: hlpr_move_subdir_files_into_target_dir.sh
#
# Description:
#   Wrapper to run move_subdir_files_into_target_dir.py inside the UV container.
#
#   The script takes a *relative path* under:
#       ${DATADIR}/1a_source/
#
#   This allows flattening:
#     - a whole dataset:
#         SC2-23
#     - or a specific subdirectory:
#         SC2-23/run_001
#
# Usage:
#   ./hlpr_move_subdir_files_into_target_dir.sh <RELATIVE_PATH_UNDER_1a_source> \
#       [--prepend-subdir] [--dry-run] [--delete-empty-subdirs]
#
# Examples:
#   ./hlpr_move_subdir_files_into_target_dir.sh SC2-24 --dry-run
#   ./hlpr_move_subdir_files_into_target_dir.sh SC2-24/Batch2_RF_possorted_genome_bam_fastq --prepend-subdir
#
# Notes:
#   - Reads DATADIR and UV_PATH from ../../config.env
#   - Passes all extra flags directly to the Python script
###############################################################################

# ----------------------------- Args -------------------------------------------

if [ $# -lt 1 ]; then
  echo "Usage: $0 <RELATIVE_PATH_UNDER_1a_source> [--prepend-subdir] [--dry-run] [--delete-empty-subdirs]"
  exit 1
fi

RELATIVE_PATH="$1"
shift  # Remaining args are passed through to Python

# ----------------------------- Resolve paths ----------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT_PATH="${SCRIPT_DIR}/../../resources/move_subdir_files_into_target_dir.py"

# ----------------------------- Config -----------------------------------------

set -a
source "$SCRIPT_DIR/../../config.env"
set +a

for var in DATADIR UV_PATH; do
  if [ -z "${!var:-}" ]; then
    echo "ERROR: $var not set in config.env"
    exit 1
  fi
done

# ----------------------------- Target dir -------------------------------------

TARGET_DIR="${DATADIR}/1a_source/${RELATIVE_PATH}"

if [ ! -d "$TARGET_DIR" ]; then
  echo "ERROR: Target directory does not exist:"
  echo "  $TARGET_DIR"
  exit 1
fi

if [ ! -f "$SCRIPT_PATH" ]; then
  echo "ERROR: Python script not found: $SCRIPT_PATH"
  exit 1
fi

# ----------------------------- Execution --------------------------------------

echo "Running move_subdir_files_into_target_dir.py"
echo "  Target directory: $TARGET_DIR"
echo "  Container:        $UV_PATH"
echo "  Extra args:       $*"
echo

singularity exec "$UV_PATH" \
  uv run python "$SCRIPT_PATH" "$TARGET_DIR" "$@"

echo
echo "Done."
