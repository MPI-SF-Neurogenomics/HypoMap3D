#!/bin/bash
set -euo pipefail

# -----------------------------------------------------------------------------
# 1ab_submit_all_samples_10x.sh
# -----------------------------------------------------------------------------
# Submits Cellranger count jobs via SLURM for all samples found in a dataset.
#
# Required config variables (from ../../config.env):
#   DATADIR
#   CELLRANGER_PATH
#   SLURM_EMAIL
#   CPU_PARTITION
#
# Positional arguments:
#   1) dataset_dir   (e.g. SC2-23)
#   2) sequencing type: sn | sc
#   3) chemistry type: auto | multiome
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
        echo "$(date): ERROR: $var not set in config.toml"
        exit 1
    fi
done

############################
# Parse arguments
############################

if [ $# -ne 3 ]; then
    echo "Usage: $0 <dataset_dir> <sn|sc> <auto|multiome>"
    exit 1
fi

DATASET_DIR=$1
SEQ_TYPE=$2
CHEMISTRY=$3

############################
# Define paths
############################

FASTQ_PATH="${DATADIR}/1b_bcl_FASTQ/${DATASET_DIR}"
OUTDIR="${DATADIR}/1c_aligned/${DATASET_DIR}"
LOGDIR="${DATADIR}/1c_aligned/logs/${DATASET_DIR}"
SLURM_LOGDIR="${LOGDIR}/slurm"

mkdir -p "$SLURM_LOGDIR"

############################
# Validate inputs
############################

if [ ! -d "$FASTQ_PATH" ]; then
    echo "ERROR: FASTQ directory not found: $FASTQ_PATH"
    exit 2
fi

############################
# Discover samples
############################

echo "Discovering samples in $FASTQ_PATH..."

SAMPLES=$(find "$FASTQ_PATH" \( -type f -o -type l \) -name "*.fastq.gz" -print0 \
    | xargs -0 -n 1 basename \
    | sed 's/_[Ss][0-9]\+_L[0-9]\+_R[0-9]\+_001.fastq.gz//' \
    | sort -u)

if [ -z "$SAMPLES" ]; then
    echo "No FASTQ files found in $FASTQ_PATH"
    exit 0
fi

echo "Found samples:"
printf "  %s\n" $SAMPLES

############################
# Submit jobs
############################

for SAMPLE in $SAMPLES; do

    RESULT_FILE="${OUTDIR}/${SAMPLE}/outs/filtered_feature_bc_matrix.h5"
    if [ -f "$RESULT_FILE" ]; then
        echo "Sample $SAMPLE already processed — skipping."
        continue
    fi

    echo "Submitting Cellranger job for sample: $SAMPLE"

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

done

echo "All eligible samples have been submitted."
