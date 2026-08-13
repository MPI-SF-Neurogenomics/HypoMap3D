#!/bin/bash

# Description:
# This script downloads Ensembl genome FASTA and GTF files,
# removes gene Gm28040 from the GTF, and logs all actions.
# Configuration is sourced from ../../config.toml.
# 
# Usage:
# ./0a_download_index_genome_files.sh
#
# Notes: 
# This script is dependent on the 'wget' utility for file download,
# and 'awk' for text processing. Ensure these utilities are available 
# on your system. The script will fail if the files to download are not 
# available on the specified URLs or if the specified Ensembl Gene ID 
# is not found in the GTF file. It will attempt to create the specified 
# directories if they do not exist. Ensure the user has write permission to the specified paths.


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

############################
# Define directories
############################

TARGET_DIR="${DATADIR}/1a_source/index_genomes/shared"
LOG_DIR="${DATADIR}/1a_source/logs"

# Ensure directories exist
mkdir -p "$TARGET_DIR"
mkdir -p "$LOG_DIR"

############################
# Start logging
############################

LOG_FILE="${LOG_DIR}/0a_download_index_genome.log"
exec > >(tee -a "$LOG_FILE")
exec 2> >(tee -a "$LOG_FILE")

echo "$(date): Starting genome download and processing"
echo "$(date): Using DATADIR = $DATADIR"
echo "$(date): Target directory = $TARGET_DIR"

############################
# URLs
############################

FILE_URL_1="https://ftp.ensembl.org/pub/release-110/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.toplevel.fa.gz"
FILE_URL_2="https://ftp.ensembl.org/pub/release-110/gtf/mus_musculus/Mus_musculus.GRCm39.110.gtf.gz"

cd "$TARGET_DIR"

############################
# Download files
############################

echo "$(date): Downloading FASTA..."
wget -q "$FILE_URL_1"

echo "$(date): Downloading GTF..."
wget -q "$FILE_URL_2"

FILE_NAME_1=$(basename "$FILE_URL_1")
FILE_NAME_2=$(basename "$FILE_URL_2")

############################
# Validate downloads
############################

if [ ! -f "$FILE_NAME_1" ]; then
    echo "$(date): ERROR: FASTA download failed"
    exit 2
fi

if [ ! -f "$FILE_NAME_2" ]; then
    echo "$(date): ERROR: GTF download failed"
    exit 3
fi

############################
# Unzip files
############################

echo "$(date): Unzipping FASTA..."
gunzip -f "$FILE_NAME_1"

echo "$(date): Unzipping GTF..."
gunzip -f "$FILE_NAME_2"
FILE_NAME_2="${FILE_NAME_2%.gz}"

############################
# Remove gene from GTF
############################

GENE_NAME="Gm28040"
OUTPUT_GTF="${FILE_NAME_2%.gtf}.woGm28040.gtf"

echo "$(date): Removing gene $GENE_NAME from GTF..."
grep -v "$GENE_NAME" "$FILE_NAME_2" > "$OUTPUT_GTF"

rm "$FILE_NAME_2"

echo "$(date): Finished successfully"