#!/bin/bash
set -euo pipefail

###############################################################################
# 1_CellBender.sh
###############################################################################
# Purpose:
#   - Process single-cell RNA-seq datasets (10x Genomics / STARsolo / Smart-seq)
#   - For 10x/STARsolo datasets:
#       * Run CellBender remove-background to remove ambient RNA
#       * Convert the resulting matrix to h5ad format for downstream analysis
#   - For Smart-seq datasets:
#       * Skip CellBender entirely
#       * Convert the raw matrix directly to h5ad (no barcode file)
#
# Arguments:
#   1) SAMPLE (string)
#       - Sample identifier in the format {dataset_id}_{sample_id}, e.g., SC2-13_SAMN12345678
#   2) Optional flag: --cap-cells
#       - Caps the estimated number of cells at 10,000 for CellBender
#
# Behavior:
#   - Detects the technology from dataset_metadata.csv
#   - Chooses the correct matrix file (Cell Ranger or STARsolo)
#   - Sets up directories for CellBender output, logs, and h5ad outputs
#   - Logs everything (stdout + stderr) to a timestamped log file in logs/<SAMPLE>/
#   - Normalizes behavior depending on technology:
#       * Smart-seq: skip CellBender, direct h5 -> h5ad conversion
#       * 10x/STARsolo: run CellBender, then h5 -> h5ad conversion
#
# Requirements:
#   - config.env must define:
#       DATADIR         : Base directory containing aligned matrices
#       CELLBENDER_PATH : Path to CellBender Singularity container
#       UV_PATH         : Path to h5 -> h5ad conversion Singularity container
#   - dataset_metadata.csv with columns: ReferenceID, Technology
#   - combined_alignment_qc_metrics_summary.csv with columns: Sample, Estimated Number of Cells
#   - h5_to_h5ad.py script in resources/ directory
#
# Logging:
#   - Logs all steps, warnings, and errors to a timestamped log file per sample
#   - Outputs also printed to stdout
#
# Example usage:
#   ./1_CellBender.sh SC2-13_SAMN12345678
#   ./1_CellBender.sh SC2-13_SAMN12345678 --cap-cells
###############################################################################

############################
# Arguments
############################

SAMPLE=$1

# Optional flag
CAP_CELLS=false
if [ "${2:-}" == "--cap-cells" ]; then
    CAP_CELLS=true
fi

############################
# Paths
############################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
set -a
source "$SCRIPT_DIR/../../config.env"
set +a

DATA_DIR="${DATADIR}/1c_aligned"
CELLBENDER_BASE="${DATADIR}/2a_cellbender"
RESOURCES_DIR="${SCRIPT_DIR}/../../resources"

DATASET_ID="${SAMPLE%%_*}"

METRICS_CSV="${DATA_DIR}/combined_alignment_qc_metrics_summary.csv"
METADATA_CSV="${SCRIPT_DIR}/../../resources/dataset_metadata.csv"

OUTDIR="${CELLBENDER_BASE}/outs/${SAMPLE}"
LOGDIR="${CELLBENDER_BASE}/logs/${SAMPLE}"
H5AD_DIR="${CELLBENDER_BASE}/outs/wo_empty_droplets"

OUTPUT_H5="${OUTDIR}/output.h5"
H5AD_PATH="${H5AD_DIR}/${SAMPLE}.h5ad"

mkdir -p "$OUTDIR" "$LOGDIR" "$H5AD_DIR"

############################
# Logging setup (everything goes to file + stdout)
############################

TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOGFILE="${LOGDIR}/${SAMPLE}_cellbender_${TIMESTAMP}.log"
exec > >(tee -a "$LOGFILE") 2>&1

echo "===================="
echo "Script started: $(date)"
echo "SAMPLE: $SAMPLE"
echo "DATASET_ID: $DATASET_ID"
echo "OUTDIR: $OUTDIR"
echo "H5AD_DIR: $H5AD_DIR"
echo "===================="

############################
# Detect matrix file
############################

MATRIX_FILE="${DATA_DIR}/${DATASET_ID}/${SAMPLE}/outs/raw_feature_bc_matrix.h5"
STARSOLO_MATRIX_FILE="${DATA_DIR}/${DATASET_ID}/${SAMPLE}/Solo.out/Gene/Unique_matrix.h5"

if [ -f "$MATRIX_FILE" ]; then
    echo "Detected Cell Ranger matrix: $MATRIX_FILE"
elif [ -f "$STARSOLO_MATRIX_FILE" ]; then
    MATRIX_FILE="$STARSOLO_MATRIX_FILE"
    echo "Detected STARsolo matrix: $MATRIX_FILE"
else
    echo "ERROR: No valid matrix file found for sample $SAMPLE"
    exit 1
fi

############################
# Determine technology
############################

IS_SMARTSEQ=false
TECHNOLOGY="Unknown"

if [ -f "$METADATA_CSV" ]; then
    TECHNOLOGY=$(awk -F',' -v dataset="$DATASET_ID" '
        NR==1 {for (i=1;i<=NF;i++){if($i=="ReferenceID") ref=i;if($i=="Technology") tech=i}}
        NR>1 && $ref==dataset {gsub(/"/,"",$tech); print $tech; exit}
    ' "$METADATA_CSV")

    if [ -n "$TECHNOLOGY" ] && echo "$TECHNOLOGY" | grep -qiE 'smart[- ]*seq'; then
        IS_SMARTSEQ=true
    fi
else
    echo "WARNING: Dataset metadata CSV not found: $METADATA_CSV"
fi

echo "Dataset technology: $TECHNOLOGY"

############################
# Smart-seq handling
############################

if $IS_SMARTSEQ; then
    echo "Smart-seq dataset detected — skipping CellBender."
    
    if [ ! -f "$MATRIX_FILE" ]; then
        echo "ERROR: Matrix file not found for Smart-seq sample: $MATRIX_FILE"
        exit 1
    fi

    if [ -f "$H5AD_PATH" ]; then
        echo "h5ad output already exists: $H5AD_PATH"
    else
        echo "Converting Smart-seq matrix to h5ad (no -b flag)..."
        cd "$(dirname "$MATRIX_FILE")" || exit 1

        singularity run --nv "$UV_PATH" \
            python -u "$RESOURCES_DIR/h5_to_h5ad.py" \
                "$MATRIX_FILE" \
                -O "$H5AD_PATH"

        echo "h5ad conversion complete: $H5AD_PATH"
    fi

    echo "Smart-seq processing done."
    echo "===================="
    exit 0
fi

############################
# Normal 10x/STARsolo handling
############################

# Extract expected cell count
RAW_CELL_COUNT=$(awk -F',' -v sample="$SAMPLE" '
    NR==1 {for (i=1;i<=NF;i++){if($i=="Sample") sample_col=i;if($i=="Estimated Number of Cells") est_col=i}}
    NR>1 && $sample_col==sample {gsub(/"/,"",$est_col); print $est_col; exit}
' "$METRICS_CSV")

if [ -z "$RAW_CELL_COUNT" ]; then
    echo "ERROR: Sample $SAMPLE not found in metrics CSV"
    exit 1
fi

if $CAP_CELLS && [ "$RAW_CELL_COUNT" -gt 10000 ]; then
    N_EXP_CELLS=10000
    echo "Estimated cells capped at 10,000"
else
    N_EXP_CELLS=$RAW_CELL_COUNT
fi

echo "Expected cells: $N_EXP_CELLS"

# Run CellBender
if [ -f "$OUTPUT_H5" ]; then
    echo "CellBender output already exists: $OUTPUT_H5"
else
    echo "Running CellBender..."
    cd "$OUTDIR"

    singularity run --nv "$CELLBENDER_PATH" \
        cellbender remove-background \
            --input "$MATRIX_FILE" \
            --output output.h5 \
            --expected-cells "$N_EXP_CELLS" \
            --model full \
            --fpr 0 \
            --cuda || { echo "CellBender failed"; exit 1; }
fi

# Convert to h5ad
if [ -f "$H5AD_PATH" ]; then
    echo "h5ad output already exists: $H5AD_PATH"
else
    echo "Converting matrix to h5ad..."
    cd "$OUTDIR"

    singularity run --nv "$UV_PATH" \
        python -u "$RESOURCES_DIR/h5_to_h5ad.py" \
            "$MATRIX_FILE" \
            -b output_cell_barcodes.csv \
            -O "$H5AD_PATH"

    echo "h5ad conversion complete: $H5AD_PATH"
fi

echo "Processing complete for sample $SAMPLE"
echo "===================="