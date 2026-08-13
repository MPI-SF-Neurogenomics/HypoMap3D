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

STAR \
  --runMode alignReads \
  --runThreadN "$THREADS" \
  --quantMode GeneCounts \
  --soloType SmartSeq \
  --soloStrand Unstranded \
  --soloUMIdedup NoDedup \
  --outFilterScoreMin 30 \
  --outSAMtype BAM Unsorted \
  --outFileNamePrefix "$OUTDIR" \
  --readFilesCommand zcat \
  --genomeDir "$TRANSCRIPTOME" \
  --sjdbGTFfile "$GTF" \
  --readFilesManifest "$MANIFEST" || exit 1

# Add arg --soloFeatures GeneFull_Ex50pAS when processing sn-RNAseq datats to include introns