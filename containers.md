# Development Containers

The use of containers provides an identical OS environment for all developers of a project - thereby avoiding the *"it works on my machine"* dilemma. Containers has consequently become gold-standard for reproducible research and smoother development. Limitations still exist however, as hardware differences can still produce differences particularly with GPU computation.

## Docker Setup

The use of docker requires admin privileges as such is highly incompatible with High Perforamce Compute (HPC) environments. As such use of `apptainer` or `singularity` is preferred.

If you wish to explore locally then please refer to docker documentation for instruction ([docs.docker.com](https://docs.docker.com/get-started/)).

## Singularity Setup (HPC Supported)

```bash
#e.g. pulling uv-r-cuda image from docker to Apptainer SIF image via singularity 
singularity pull uv-r-cuda_0.2.1.sif docker://ehvr20/uv-r-cuda:0.2.1

# Running the image interactively
singularity exec uv-r-cuda_0.2.1.sif bash

# Running `uv` within the image
singularity exec uv-r-cuda_0.2.1.sif uv --help
```

## UV-R-CUDA (ehvr20/uv-r-cuda:0.2.1)

This image is a highly flexible, lightweight image comprising both `uv` and `R` alongside required `cuda` dependancies. Python environments are highly reproducible thanks to `uv` as such we use the container environment to provide a stable OS and wheel build behaviour across all developers.

## STARsolo (ehvr20/star:2.7.11b)

A lightweight image for running STAR and STARsolo (2.7.11b). Used for alignment of raw 10X technology reads from FASTQs into count matrices (.h5/mtx).

## Cellranger (ehvr20/cellranger:9.0.0)

A lightweight image for running 10X Cellranger (9.0.0). Used for alignment of raw drop-seq and smart-seq technology reads from FASTQs into count matrices (.h5/mtx).

## Cellbender (us.gcr.io/broad-dsde-methods/cellbender:0.3.0)

The offical cellbender docker image (v0.3.0) for more details please refer to [gh/broadinstitute/CellBender](https://github.com/broadinstitute/CellBender) 

## sra-tools (ncbi/sra-tools)

The offical NCBI SRA toolkit docker image. For more details please refer to [gh/ncbi/sra-tools](https://github.com/ncbi/sra-tools/)

## r_scvi_docker_v4 (lsteuernagel/r_scvi_docker_v4:v9)
A fixed R environment providing a stable sourcing of required R packages such as Seurat.