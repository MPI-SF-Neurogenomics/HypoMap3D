#!/bin/bash
set -euo pipefail

###############################################################################
# 1a_submit_all_samples_cellbender.sh
###############################################################################
# Submits SLURM jobs to process all single-cell datasets:
#   - Runs CellBender for 10x/STARsolo datasets
#   - Skips CellBender for Smart-seq datasets, but converts to h5ad
#
# Sample discovery:
#   - All samples under ${DATADIR}/1c_aligned/<dataset>/<sample>/
#   - Sample format: {dataset_id}_{sample_id}
#
# Skipping logic:
#   - If ${DATADIR}/2a_cellbender/outs/<sample>/output.h5 exists AND h5ad exists,
#     the sample is skipped
#
# Optional argument:
#   --cap-cells : caps estimated cells at 10,000 for CellBender
###############################################################################

############################
# Arguments
############################
CAP_CELLS_FLAG=${1:-}   # Optional: --cap-cells

############################
# Resolve paths
############################
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

############################
# Load configuration
############################
set -a
source "$SCRIPT_DIR/../../config.env"
set +a

for var in DATADIR CELLBENDER_PATH SLURM_EMAIL GPU_PARTITION; do
    if [ -z "${!var:-}" ]; then
        echo "$(date): ERROR: $var not set in config.env"
        exit 1
    fi
done

############################
# Fixed directories
############################
DATA_DIR="${DATADIR}/1c_aligned"
CELLBENDER_BASE="${DATADIR}/2a_cellbender"
SLURMLOG_BASE="${CELLBENDER_BASE}/logs"

mkdir -p "$SLURMLOG_BASE"

############################
# Discover samples
############################
echo "Discovering samples in: $DATA_DIR"

SAMPLES=()

for DATASET_DIR in "$DATA_DIR"/*; do
    [ -d "$DATASET_DIR" ] || continue

    for SAMPLE_DIR in "$DATASET_DIR"/*; do
        [ -d "$SAMPLE_DIR" ] || continue
        SAMPLE=$(basename "$SAMPLE_DIR")
        SAMPLES+=("$SAMPLE")
    done
done

if [ "${#SAMPLES[@]}" -eq 0 ]; then
    echo "No samples found. Exiting."
    exit 0
fi

############################
# Submit jobs
############################
for SAMPLE in "${SAMPLES[@]}"; do

    OUTPUT_H5="${CELLBENDER_BASE}/outs/${SAMPLE}/output.h5"
    H5AD_PATH="${CELLBENDER_BASE}/outs/wo_empty_droplets/${SAMPLE}.h5ad"

    # Skip only if BOTH CellBender output and h5ad already exist
    if [ -f "$OUTPUT_H5" ] && [ -f "$H5AD_PATH" ]; then
        echo "Skipping $SAMPLE (CellBender + h5ad already exist)"
        continue
    fi

    echo "Submitting processing job for sample: $SAMPLE"

    SAMPLE_LOGDIR="${SLURMLOG_BASE}/${SAMPLE}/slurm"
    mkdir -p "$SAMPLE_LOGDIR"

    sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=cellbender_${SAMPLE}
#SBATCH --output=${SAMPLE_LOGDIR}/${SAMPLE}_slurm.out
#SBATCH --error=${SAMPLE_LOGDIR}/${SAMPLE}_slurm.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=72G
#SBATCH --partition=${GPU_PARTITION}
#SBATCH --mail-type=END,FAIL,CANCEL
#SBATCH --mail-user=${SLURM_EMAIL}
#SBATCH --nice

set -euo pipefail

# Run the main processing script; it handles Smart-seq and 10x datasets
srun "$SCRIPT_DIR/1_CellBender.sh" "$SAMPLE" ${CAP_CELLS_FLAG}
EOF

done

echo "All eligible jobs submitted."