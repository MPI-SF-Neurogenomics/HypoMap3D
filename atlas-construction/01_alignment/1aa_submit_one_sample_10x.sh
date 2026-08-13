#!/bin/bash
set -euo pipefail

# -----------------------------------------------------------------------------
# 1aa_submit_one_sample_10x.sh
# -----------------------------------------------------------------------------
# Submits a single Cellranger count job via SLURM for one sample.
#
# Required config variables (from ../../config.env):
#   DATADIR
#   CELLRANGER_PATH
#   SLURM_EMAIL
#   CPU_PARTITION
#
# Positional arguments:
#   1) dataset_dir   (e.g. SC2-23)
#   2) sample name   (e.g. SC2-24_SAMN43828891)
#   3) sequencing type: sn | sc
#   4) chemistry type: auto | multiome
# -----------------------------------------------------------------------------

############################
# Load configuration
############################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

set -a
source "$SCRIPT_DIR/../../config.env"
set +a

for var in DATADIR CELLRANGER_PATH SLURM_EMAIL CPU_PARTITION; do
    if [ -z "${!var:-}" ]; then
        echo "$(date): ERROR: $var not set in config.env"
        exit 1
    fi
done

############################
# Parse arguments
############################

if [ $# -ne 4 ]; then
    echo "Usage: $0 <dataset_dir> <sample> <sn|sc> <auto|multiome>"
    exit 1
fi

DATASET_DIR=$1
SAMPLE=$2
SEQ_TYPE=$3
CHEMISTRY=$4

############################
# Define paths
############################

FASTQ_PATH="${DATADIR}/1b_bcl_FASTQ/${DATASET_DIR}"
OUTDIR="${DATADIR}/1c_aligned/${DATASET_DIR}"
LOGDIR="${DATADIR}/1c_aligned/logs/${DATASET_DIR}"
SLURM_LOGDIR="${LOGDIR}/slurm"

mkdir -p "$SLURM_LOGDIR"

############################
# Skip if already completed
############################

RESULT_FILE="${OUTDIR}/${SAMPLE}/outs/filtered_feature_bc_matrix.h5"
if [ -f "$RESULT_FILE" ]; then
    echo "Cellranger has already been run on sample $SAMPLE. Skipping submission."
    exit 0
fi

############################
# Validate inputs
############################

if [ ! -d "$FASTQ_PATH" ]; then
    echo "ERROR: FASTQ directory not found: $FASTQ_PATH"
    exit 2
fi

if [ ! -f "$CELLRANGER_PATH" ]; then
    echo "ERROR: Cellranger Singularity image not found: $CELLRANGER_PATH"
    exit 3
fi

############################
# Submit SLURM job
############################

echo "Submitting Cellranger job for sample: $SAMPLE"
echo "Dataset: $DATASET_DIR"

sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=cellranger_${SAMPLE}
#SBATCH --output=${SLURM_LOGDIR}/${SAMPLE}_slurm.out
#SBATCH --error=${SLURM_LOGDIR}/${SAMPLE}_slurm.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=48G
#SBATCH --partition=${CPU_PARTITION}
#SBATCH --mail-type=END,FAIL,CANCEL
#SBATCH --mail-user=${SLURM_EMAIL}

set -euo pipefail

srun "${SCRIPT_DIR}/1a_10x.sh" \
    "${DATASET_DIR}" \
    "${SAMPLE}" \
    "${SEQ_TYPE}" \
    "${CHEMISTRY}"
EOF
