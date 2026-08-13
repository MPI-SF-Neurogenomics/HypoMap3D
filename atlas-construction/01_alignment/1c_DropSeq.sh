#!/bin/bash

set -euo pipefail

OUTDIR='results/'
THREADS=56

# Indexed transcriptome is generated using scripts in 'Reference_Genome'
TRANSCRIPTOME='path_to_transcriptome_dir'
GTF='path_to_gtf_file'

# A tsv file of format: Read2\tRead1\tSampleID
# N/B If no Read1 replace with a '-' instead
MANIFEST='path_to_fastq_manifest_file'

# Cell Barcode (CB) and UMI lengths must be manually defined.
CB_LEN=12
CB_START=0
UMI_LEN=8
UMISTART=13

STAR \
 --runMode alignReads \
 --runThreadN "$THREADS" \
 --quantMode GeneCounts \
 --soloType CB_UMI_Simple \
 --soloStrand Forward \
 --soloCBstart "$CB_START" --soloCBlen "$CB_LEN" --soloUMIstart "$UMI_START" --soloUMIlen "$UMI_LEN" --soloCBwhitelist None --soloBarcodeReadLength 0 \
 --soloCellFilter EmptyDrops_CR \
 --soloUMIdedup 1MM_CR \
 --limitOutSJcollapsed 20000000 \
 --outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3 \
 --outSAMtype BAM Unsorted \
 --outFileNamePrefix "$OUTDIR" \
 --readFilesCommand zcat \
 --genomeDir "$TRANSCRIPTOME" \
 --sjdbGTFfile "$GTF" \
 --readFilesManifest "$MANIFEST"

# Add arg --soloFeatures GeneFull_Ex50pAS when processing sn-RNAseq datats to include introns