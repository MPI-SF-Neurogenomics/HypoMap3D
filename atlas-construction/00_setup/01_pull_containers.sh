#!/bin/bash
set -euo pipefail


# Check if a path was provided
if [ -z "$1" ]; then
  echo "Usage: $0 /path/to/data/root/directory"
  exit 1
fi

DATA_ROOT_DIR="$1"

singularity pull "${DATA_ROOT_DIR}/uv_r_cuda_0.2.1.sif docker://ehvr20/uv-r-cuda:0.2.1"
singularity pull "${DATA_ROOT_DIR}/STAR_2.7.11b.sif docker://ehvr20/star:2.7.11b"
singularity pull "${DATA_ROOT_DIR}/cellranger_9.0.0.sif docker://ehvr20/cellranger:9.0.0"
singularity pull "${DATA_ROOT_DIR}/cellbender_0.3.0.sif docker://us.gcr.io/broad-dsde-methods/cellbender:0.3.0"
singularity pull "${DATA_ROOT_DIR}/r_scvi_docker_v4_v9 docker://lsteuernagel/r_scvi_docker_v4:v9"