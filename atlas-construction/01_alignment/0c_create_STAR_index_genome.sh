#!/bin/bash
set -euo pipefail

# -----------------------------------------------------------------------------
# 0c_misc_STAR_index_genome.sh
# -----------------------------------------------------------------------------
# Generates a STAR genome index for GRCm39 using a Singularity container.
# Configuration is sourced from ../../config.toml.
#
# Required config variables:
#   DATADIR
#   STAR_PATH
# -----------------------------------------------------------------------------

############################
# Load configuration
############################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"set -a

source "$SCRIPT_DIR/../../config.env"
set +a

if [ -z "${DATADIR:-}" ]; then
    echo "$(date): ERROR: DATADIR not set in config.env"
    exit 1
fi

if [ -z "${STAR_PATH:-}" ]; then
    echo "$(date): ERROR: STAR_PATH not set in config.env"
    exit 1
fi

############################
# Define directories
############################

FILES_DIR="${DATADIR}/1a_source/index_genomes/shared"
GENOME_BASE_DIR="${DATADIR}/1a_source/index_genomes/star_index_genome"
GENOME_DIR="${GENOME_BASE_DIR}/GRCm39/star-v2.7.11b"
LOG_DIR="${DATADIR}/1a_source/logs"

# Ensure required directories exist
mkdir -p "$GENOME_DIR"
mkdir -p "$LOG_DIR"

############################
# Start logging
############################

LOG_FILE="${LOG_DIR}/0c_misc_STAR_index_genome.log"
exec > >(tee -a "$LOG_FILE")
exec 2> >(tee -a "$LOG_FILE")

echo "$(date): Starting STAR genome indexing"
echo "$(date): DATADIR = $DATADIR"
echo "$(date): Files directory = $FILES_DIR"
echo "$(date): Genome directory = $GENOME_DIR"
echo "$(date): Singularity image = $STAR_PATH"

############################
# Define input files
############################

FASTA="${FILES_DIR}/Mus_musculus.GRCm39.dna.toplevel.fa"
GTF="${FILES_DIR}/Mus_musculus.GRCm39.110.woGm28040.gtf"

# Validate inputs
if [ ! -f "$FASTA" ]; then
    echo "$(date): ERROR: FASTA file not found: $FASTA"
    exit 2
fi

if [ ! -f "$GTF" ]; then
    echo "$(date): ERROR: GTF file not found: $GTF"
    exit 3
fi

if [ ! -f "$STAR_PATH" ]; then
    echo "$(date): ERROR: Singularity image not found: $STAR_PATH"
    exit 4
fi

############################
# Run STAR genomeGenerate
############################

echo "$(date): Running STAR genomeGenerate..."

singularity exec "$STAR_PATH" STAR \
    --runThreadN 56 \
    --runMode genomeGenerate \
    --genomeDir "$GENOME_DIR" \
    --genomeFastaFiles "$FASTA" \
    --sjdbGTFfile "$GTF"

echo "$(date): STAR genome index generation completed successfully"
