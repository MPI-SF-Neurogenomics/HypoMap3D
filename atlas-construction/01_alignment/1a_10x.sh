#!/bin/bash
set -euo pipefail

# -----------------------------------------------------------------------------
# 1a_10x.sh
# -----------------------------------------------------------------------------
# Runs cellranger count using Singularity.
#
# Required config variables (from ../../config.env):
#   DATADIR
#   CELLRANGER_PATH
#
# Positional arguments:
#   1) dataset_dir   (e.g. SC2-23)
#   2) sample name
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

if [ -z "${DATADIR:-}" ]; then
    echo "$(date): ERROR: DATADIR not set in config.env"
    exit 1
fi

if [ -z "${CELLRANGER_PATH:-}" ]; then
    echo "$(date): ERROR: CELLRANGER_PATH not set in config.env"
    exit 1
fi

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

FASTQ_PATH="${DATADIR}/1c_bcl_FASTQ/${DATASET_DIR}"
OUTDIR="${DATADIR}/1c_aligned/${DATASET_DIR}"
LOGDIR="${DATADIR}/1c_aligned/logs/${DATASET_DIR}"
TRANSCRIPTOME="${DATADIR}/1a_source/index_genomes/cellranger_index_genome"

mkdir -p "$OUTDIR" "$LOGDIR"

TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOGFILE="${LOGDIR}/${SAMPLE}_cellranger_run_${TIMESTAMP}.log"

############################
# Start logging
############################

exec > >(tee -a "$LOGFILE") 2>&1

echo "Logging started at $(date)"
echo "DATADIR: $DATADIR"
echo "CELLRANGER_IMAGE: $CELLRANGER_PATH"
echo "DATASET_DIR: $DATASET_DIR"
echo "FASTQ_PATH: $FASTQ_PATH"
echo "OUTDIR: $OUTDIR"
echo "LOGDIR: $LOGDIR"
echo "SAMPLE: $SAMPLE"
echo "SEQ_TYPE: $SEQ_TYPE"
echo "CHEMISTRY: $CHEMISTRY"

############################
# Validate inputs
############################

if [ ! -d "$FASTQ_PATH" ]; then
    echo "ERROR: FASTQ directory not found: $FASTQ_PATH"
    exit 2
fi

if [ ! -d "$TRANSCRIPTOME" ]; then
    echo "ERROR: Transcriptome directory not found: $TRANSCRIPTOME"
    exit 3
fi

if [ ! -f "$CELLRANGER_PATH" ]; then
    echo "ERROR: Cellranger Singularity image not found: $CELLRANGER_PATH"
    exit 4
fi

############################
# Sequencing mode logic
############################

case "$SEQ_TYPE" in
    sn) INCLUDE_INTRONS=true ;;
    sc) INCLUDE_INTRONS=false ;;
    *)
        echo "ERROR: SEQ_TYPE must be 'sn' or 'sc'"
        exit 5
        ;;
esac

case "$CHEMISTRY" in
    auto) CHEMISTRY_OPTION="auto" ;;
    multiome) CHEMISTRY_OPTION="ARC-v1" ;;
    *)
        echo "ERROR: CHEMISTRY must be 'auto' or 'multiome'"
        exit 6
        ;;
esac

############################
# Skip if already completed
############################

RESULT_FILE="${OUTDIR}/${SAMPLE}/outs/filtered_feature_bc_matrix.h5"
if [ -f "$RESULT_FILE" ]; then
    echo "Cellranger output already exists for sample $SAMPLE — skipping."
    exit 0
fi

############################
# Run Cellranger
############################

cd "$OUTDIR"

echo "Running cellranger count..."

singularity run "$CELLRANGER_PATH" \
    cellranger count \
    --transcriptome "$TRANSCRIPTOME" \
    --fastqs "$FASTQ_PATH" \
    --sample "$SAMPLE" \
    --id "$SAMPLE" \
    --include-introns "$INCLUDE_INTRONS" \
    --chemistry "$CHEMISTRY_OPTION"

echo "Cellranger count completed for sample $SAMPLE"
