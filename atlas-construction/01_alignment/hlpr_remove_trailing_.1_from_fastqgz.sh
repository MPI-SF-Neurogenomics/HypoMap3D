#!/bin/bash
set -euo pipefail

###############################################################################
# Script Name: hlpr_remove_trailing_.1_from_fastqgz.sh
#
# Description:
#   This script renames all files ending with ".fastq.gz.1" by removing
#   the trailing ".1". The target directory is derived from:
#
#       ${DATADIR}/1a_source/${DATASET_DIR}
#
# Usage:
#   ./hlpr_remove_trailing_.1_from_fastqgz.sh <DATASET_DIR> [--dry-run|-n]
#
# Example:
#   ./hlpr_remove_trailing_.1_from_fastqgz.sh SC2-23
#   ./hlpr_remove_trailing_.1_from_fastqgz.sh SC2-23 --dry-run
#
# Options:
#   --dry-run | -n   Show what would be renamed without making any changes.
#
# Notes:
#   - Non-recursive: only files directly inside the dataset directory
#   - Requires DATADIR to be defined in ../../config.env
###############################################################################

# ----------------------------- Args -------------------------------------------

if [ $# -lt 1 ]; then
  echo "Usage: $0 <DATASET_DIR> [--dry-run|-n]"
  exit 1
fi

DATASET_DIR="$1"
DRY_RUN=false

if [[ "${2:-}" == "--dry-run" || "${2:-}" == "-n" ]]; then
  DRY_RUN=true
fi

# ----------------------------- Config -----------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

set -a
source "$SCRIPT_DIR/../../config.env"
set +a

if [ -z "$DATADIR" ]; then
  echo "ERROR: DATADIR not set in config.toml"
  exit 1
fi

# ----------------------------- Paths ------------------------------------------

TARGET_DIR="${DATADIR}/1a_source/${DATASET_DIR}"

if [ ! -d "$TARGET_DIR" ]; then
  echo "ERROR: '$TARGET_DIR' is not a valid directory."
  exit 1
fi

# ----------------------------- Execution --------------------------------------

echo "Processing directory: $TARGET_DIR"
echo "Dry run: $DRY_RUN"

shopt -s nullglob

for file in "$TARGET_DIR"/*.fastq.gz.1; do
  new_name="${file%.1}"

  if [ "$DRY_RUN" = true ]; then
    echo "[DRY RUN] Would rename: $(basename "$file") -> $(basename "$new_name")"
  else
    mv "$file" "$new_name"
    echo "Renamed: $(basename "$file") -> $(basename "$new_name")"
  fi
done

echo "Done."
