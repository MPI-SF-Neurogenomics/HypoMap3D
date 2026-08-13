#!/bin/bash

# : '
#     Submits a SLURM job to run the final marker gene detection pipeline.
#
#     Usage:
#         ./6_submit_final_markers.sh 
# '

set -euo pipefail

# ----------------------------- Paths -----------------------------
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
PYTHON_SCRIPT="$SCRIPT_DIR/5_run_final_markers.py"

# ----------------------------- Config -----------------------------

set -a
source "$SCRIPT_DIR/../../config.env"
set +a

LOGDIR="$DATADIR/5_clustering/logs"
mkdir -p "$LOGDIR"

TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
JOB_NAME="final_markers"
LOGFILE="$LOGDIR/${JOB_NAME}_${TIMESTAMP}"


# ----------------------------- Command Build -----------------------------
SBATCH_COMMAND="uv run $PYTHON_SCRIPT"

# ----------------------------- Submit Job -----------------------------
sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=${JOB_NAME}
#SBATCH --output=${LOGFILE}_slurm.out
#SBATCH --error=${LOGFILE}_slurm.err
#SBATCH --time=300:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=100G
#SBATCH --partition=${CPU_PARTITION}
#SBATCH --nice
#SBATCH --mail-type=END,FAIL,CANCEL
#SBATCH --mail-user=${SLURM_EMAIL}

cd "$SCRIPT_DIR"

echo "=== Final Marker Gene Analysis Execution ==="
echo "Date: \$(date)"
echo "Working Directory: \$(pwd)"
echo "IMAGE: $UV_PATH"
echo "SCRIPT: $PYTHON_SCRIPT"
echo "LOGFILE: $LOGFILE"
echo ""

singularity exec "$UV_PATH" $SBATCH_COMMAND
EOF
