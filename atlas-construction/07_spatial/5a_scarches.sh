#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

set -a
source "$SCRIPT_DIR/../../config.env"
set +a

VENV_PY="${SCRIPT_DIR}/../../.venv/bin/python3"

MODEL="${DATADIR}/5_clustering/all_control_posthoc_qc/models/scvi"

# Use an array of datasets
DATASETS=(
    "${DATADIR}/7_spatial/bounded"/SP*.h5ad
    "${DATADIR}/7_spatial/harmonized/SP2-03_no_qc.h5ad" # Open-ST dataset needs an unbounded and unwrapped variant as to create all figures
)

OUTDIR="${DATADIR}/7_spatial/embeddings"
mkdir -p "$OUTDIR"

LOGDIR="$OUTDIR/logs"
mkdir -p "$LOGDIR"

for DATASET in "${DATASETS[@]}"; do
    BASENAME="$(basename "$DATASET" .h5ad)"

    # special naming for the unbounded dataset
    if [[ "$DATASET" == *"/harmonized/"* ]]; then
        OUTFILE="${OUTDIR}/${BASENAME}_unbounded_X_embed.csv.gz"
    else
        OUTFILE="${OUTDIR}/${BASENAME}_X_embed.csv.gz"
    fi

    if [[ -f "$OUTFILE" ]]; then
        echo "Skipping ${BASENAME} (output exists)"
        continue
    fi

    # ----------------------------- Metadata -----------------------------
    TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
    JOB_NAME="scArches_${BASENAME}"
    LOGFILE="$LOGDIR/${JOB_NAME}_${TIMESTAMP}"

# ----------------------------- Submit Job -----------------------------
sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=${JOB_NAME}
#SBATCH --output=${LOGFILE}_slurm.out
#SBATCH --error=${LOGFILE}_slurm.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=240G
#SBATCH --partition=${GPU_PARTITION}
#SBATCH --nice
#SBATCH --mail-type=END,FAIL,CANCEL
#SBATCH --mail-user=${SLURM_EMAIL}

cd "$SCRIPT_DIR"

echo "=== HypoMap Label Transfer Script Execution ==="
echo "Date: \$(date)"
echo "Working Directory: \$(pwd)"
echo "Dataset: $BASENAME"
echo ""

singularity exec --nv "${UV_PATH}" \
    "${VENV_PY}" -m hm3d_utils.cli.scarches \
    "$MODEL" \
    "$DATASET" \
    --outfile "$OUTFILE" \
    --epochs 24 \
    --batch_size 1024 \
    --seed 0 \
    --verbose
EOF


done