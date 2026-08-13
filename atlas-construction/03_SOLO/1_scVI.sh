#!/bin/bash

# ----------------------------- Paths -----------------------------
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

SCVI_SCRIPT="$SCRIPT_DIR/../../resources/scVI.py"
PARAMS="$SCRIPT_DIR/../../resources/parameters_solo.json"

# ----------------------------- SLURM Config -----------------------------
set -a
source "$SCRIPT_DIR/../../config.env"
set +a

LOGDIR="$DATADIR/3_solo/logs"
mkdir -p "$LOGDIR"

# Timestamped log filename
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOGFILE="$LOGDIR/scVI_${TIMESTAMP}.log"

# Redirect stdout and stderr to log file
exec > "$LOGFILE"
exec 2>&1

# cd to SCRIPT_DIR 
cd SCRIPT_DIR

# Log the configuration
echo "=== scVI Script Execution ==="
echo "Date: $(date)"
echo "Working Directory: $(pwd)"
echo "IMAGE: $UV_PATH"
echo "SCRIPT: $SCVI_SCRIPT"
echo "INPUT: $PARAMS"
echo "LOGFILE: $LOGFILE"
echo "============================="
echo ""

# Run the script inside the container
singularity exec --nv "$UV_PATH" uv run "$SCVI_SCRIPT" -v "$PARAMS"
