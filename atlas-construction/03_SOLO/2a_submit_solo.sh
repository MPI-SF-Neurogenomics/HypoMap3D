#!/bin/bash

# : '
#     This script submits a SLURM job to run SOLO doublet detection using a Singularity container.
#     To speed things up, one job is submitted for either control or non-control samples as to have a bit of parallelization.
#
#     Requirements:
#     - Assumes 2_run_solo.py script is located in the same directory as this script.
#     - SLURM logs are saved to: logs/SOLO_prediction_${MODE}_${TIMESTAMP}.log"
#
#     Example usage:
#     ./2a_submit_solo.sh 
# '

set -euo pipefail

# ----------------------------- Paths -----------------------------

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# ----------------------------- Config -----------------------------
set -a
source "$SCRIPT_DIR/../../config.env"
set +a

SOLO_SCRIPT="$SCRIPT_DIR/2_run_solo.py"
LOGDIR="$DATADIR/3_solo/logs"
mkdir -p "$LOGDIR"

# ----------------------------- cd to SCRIPT_DIR for uv venv -----------------------------

cd $SCRIPT_DIR

# ----------------------------- SLURM Job Submission -----------------------------

MODES=(control non_control)

for MODE in "${MODES[@]}"; do
    TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
    LOGFILE="$LOGDIR/solo_${MODE}_${TIMESTAMP}.log"

    echo "Submitting SOLO job for mode: $MODE"

    sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=SOLO_${MODE}
#SBATCH --output=${LOGFILE}
#SBATCH --error=${LOGFILE}
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=128G
#SBATCH --partition=${GPU_PARTITION}
#SBATCH --nice
#SBATCH --mail-type=END,FAIL,CANCEL
#SBATCH --mail-user=${SLURM_EMAIL}

cd "$SCRIPT_DIR"

echo "=== SOLO Script Execution ==="
echo "Date: \$(date)"
echo "Working Directory: \$(pwd)"
echo "IMAGE: $UV_PATH"
echo "SCRIPT: $SOLO_SCRIPT"
echo "MODE: $MODE"
echo "LOGFILE: $LOGFILE"

echo "============================="
echo ""

singularity exec --nv "$UV_PATH" uv run "$SOLO_SCRIPT" --mode "$MODE"
EOF

done