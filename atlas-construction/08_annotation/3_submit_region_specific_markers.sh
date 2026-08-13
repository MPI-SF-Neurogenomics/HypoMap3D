#!/bin/bash
# : '
# Submit region-specific marker jobs.
# One SLURM job is submitted for every combination of:
# - ANALYSIS_LEVEL in LEVELS
# - ACRO in ACROS
#
# The Python script is executed inside the UV Singularity container
# specified in config.env.
#
# Example usage:
# ./3_submit_region_specific_markers.sh
# '

set -euo pipefail

# ----------------------------- Paths -----------------------------

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# ----------------------------- Config -----------------------------

set -a
source "$SCRIPT_DIR/../../config.env"
set +a

PY_SCRIPT="$SCRIPT_DIR/3_region_specific_markers.py"

LOGDIR="$DATADIR/8_annotation/logs"
mkdir -p "$LOGDIR"

# ----------------------------- User Parameters -----------------------------

INCLUDE_NN=false

LEVELS=(
2
3
4
5
)

# all the regions that are top regions for some L3 cluster
ACROS=(
MPO
LPO
ZI
PH
LHA
PVHd
ARH
DMH
PS
MPN
PVpo
AHN
TU
SCH
VMH
TMv
PMv
SBPV
PVi
AVPV
MEPO
PVH
LM
PSTN
STN
MEA
BST
LSr
SI
CP
VTA
MMl
VAL
VPM
PMd
MMme
MS
LSv
MB
# )

# the remaining HY regions
# ACROS=(
PVa
ADP
PVp
RCH
AVP
ASO
PST
SO
SUM
VLPO
OV
FF
PD
TMd
ME
VMPO
PeF
MMm
MMp
MMd
)

# ----------------------------- Job Submission -----------------------------
for LEVEL in "${LEVELS[@]}"; do
    for ACRO in "${ACROS[@]}"; do

        TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
        LOGFILE="$LOGDIR/markers_L${LEVEL}_${ACRO}_${TIMESTAMP}.log"

        echo "Submitting: L${LEVEL} ${ACRO}"

        if [ "$INCLUDE_NN" = true ]; then
            EXTRA_ARGS="--include-nn"
            JOB_SUFFIX="_NN"
        else
            EXTRA_ARGS=""
            JOB_SUFFIX=""
        fi

        sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=mk_L${LEVEL}_${ACRO}${JOB_SUFFIX}
#SBATCH --output=${LOGFILE}
#SBATCH --error=${LOGFILE}
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=64G
#SBATCH --partition=${CPU_PARTITION}
#SBATCH --nice
#SBATCH --mail-type=FAIL,CANCEL,END
#SBATCH --mail-user=${SLURM_EMAIL}

cd "$SCRIPT_DIR"

echo "=== Region-specific markers ==="
echo "Date: \$(date)"
echo "Level: ${LEVEL}"
echo "ACRO: ${ACRO}"
echo "INCLUDE_NN: ${INCLUDE_NN}"
echo "Script: ${PY_SCRIPT}"
echo "Logfile: ${LOGFILE}"
echo ""

singularity exec "$UV_PATH" \
    uv run "$PY_SCRIPT" \
    --analysis-level ${LEVEL} \
    --acro ${ACRO} \
    ${EXTRA_ARGS}
EOF

        sleep 0.1

    done
done