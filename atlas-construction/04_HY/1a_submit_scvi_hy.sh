#!/bin/bash

# : '
#     This script submits a SLURM job to run the scVI pipeline using a Singularity container.
#
#     Requirements:
#     - Assumes scVI pipeline is triggered via ../../resources/scVI.py
#     - Assumes parameter file named parameters_HY.json exists in ../../resoures/
#     - SLURM logs are saved to: logs/scVI_<TIMESTAMP>.log
#
#     Example usage:
#     ./1a_submit_scvi_HY.sh
# '

set -euo pipefail

# ----------------------------- Paths -----------------------------
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# ----------------------------- SLURM Config -----------------------------
set -a
source "$SCRIPT_DIR/../../config.env"
set +a

SCVI_SCRIPT="$SCRIPT_DIR/../../resources/scVI.py"
PARAMS="$SCRIPT_DIR/../../resources/parameters_HY.json"
LOGDIR="$DATADIR/4_HY/logs"
mkdir -p "$LOGDIR"

# Timestamped log filename
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOGFILE="$LOGDIR/scVI_${TIMESTAMP}.log"

# ----------------------------- SLURM Job Submission -----------------------------

echo "Submitting scVI HY job"

sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=scVI_HY
#SBATCH --output=${LOGFILE}
#SBATCH --error=${LOGFILE}
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=128G
#SBATCH --partition=${GPU_PARTITION}
#SBATCH --nice
#SBATCH --mail-type=END,FAIL,CANCEL
#SBATCH --mail-user=${SLURM_EMAIL}

cd "$SCRIPT_DIR"

echo "=== scVI Script Execution ==="
echo "Date: \$(date)"
echo "Working Directory: \$(pwd)"
echo "IMAGE: $UV_PATH"
echo "SCRIPT: $SCVI_SCRIPT"
echo "INPUT: $PARAMS"
echo "LOGFILE: $LOGFILE"

echo "============================="
echo ""

singularity exec --nv "$UV_PATH" uv run "$SCVI_SCRIPT" -v "$PARAMS"
EOF
