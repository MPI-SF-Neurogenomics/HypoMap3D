# HypoMap3D

>This repository contains code for the generation and analysis of the HypoMap3D molecular atlas.

The HypoMap3D molecular atlas is an update on the original [HypoMap](https://www.nature.com/articles/s42255-022-00657-y) boasting increased granularity with the addition of further datasets and the identification of rare populations through improved clustering methodology. The integrated dataset contains >1.5 million cells from 43 datasets across the mouse hypothalamus of which 800k are considered *control* or *steady-state*. The integration of this single cell atlas with hypothalamic spatial transcriptomics data has allowed greater investigation into the spatial heterogeneity of the hypothalamus. In combination we identify 716 cell-types within the mouse hypothalamus, alongside their spatial localization.

> "HypoMap3D: An enhanced single-cell and spatial atlas of the mouse hypothalamus".  
> – available at XXX

## Data Access

* For information on how to access HypoMap3D data objects, please refer to the [HypoMap3D Documentation](https://mpi-sf-neurogenomics.github.io/HypoMap3D) hosted via GitPages.

* For the interactive exploration of the HypoMap3D please check out the HypoMap3D Portal available at XXX.

## Directory Information

- `atlas-construction`: source code for the construction of the HypoMap3D atlas
- `dev-analysis`: analyses run throughout development cycle of HypoMap3D
- `docs`: html-based documentation generated via sphinx and hosted via **GitPages**
- `example-notebooks`: notebooks running though standard analysis workflows on HypoMap3D dataset
- `publication`: source code for the generation of data, figures and tables for the HypoMap3D publication [REF]
- `resources`: general resource folder
- `src`: source code for documentation generation and HypoMap3D python packages

## Packages

- `hm3d_utils`: utility functions specific to the HypoMap3D project.
- `project_utils`: lightweight project management utilities.

## Environment

The package dependancies required to run the example analysis workflows (see [example-notebooks](./example-notebooks/)) and utilize functions defined within the `hm3d_utils` package are provided in [pyproject.toml](./pyproject.toml)

**uv installation: (highly recommended)**
```bash
uv sync
```
**pip installation:**
```bash
pip install -r requirements.txt
```

*If you wish to reconstruct the exact environment used for the development of HypoMap3D. Please refer to "[Restoring development python environment](#restoring-development-python-environment)" below.*

## HypoMap3D Development

> All scripts were developed to run exclusively on Linux operating systems. We do not guarantee compatability with Windows or macOS.

### Path Configuration

All scripts operate relative to paths defined in [config.env](./config.env). Before starting ensure these paths are correct and that the requisite singularity images have been provided (see below).

### Containerisation

All development of HypoMap3D must be carried out within the [uv-r-cuda:0.2.1](https://hub.docker.com/layers/ehvr20/uv-r-cuda/0.2.1/images/sha256-65aaf2d54fbade68c8115b7b7a96a00bdc6ca93807a3a93e14ed59cbf6743115) container or an alternative container if specified (e.g. STARsolo container). 

For more information on setting up containers as well as details on other containers used, please refer to [containers.md](./containers.md).

Scripts for pulling docker containers to singularity can also be found at [./atlas-construction/00_setup](./atlas-construction/00_setup/01_pull_containers.sh)

### Restoring development python environment

The exact package specification is detailed in the [uv.lock](./uv.lock).

To restore this exact environment, run the following:
```bash
cd HypoMap3D
uv sync --frozen
```

### Running R scripts

For running R scripts, please use the R environment provided within the `r_scvi_docker_v4:v9` docker image. For more details, please refer to [r_scvi_docker_v4](./containers.md/#r_scvi_docker_v4-lsteuernagelr_scvi_docker_v4v9).

### Generating documentation

An optional dependancy group, `docs`, is included in HypoMap3D usage environment for the generation of `sphinx` documentation in [docs/](./docs/).

Whilst running the `uv-r-cuda` container, use the following:
```bash
uv sync --frozen --group docs # restore venv
. .venv/bin/activate # enter venv
sphinx build --write-all src/_sphinx docs/ # build docs
```

## Notes

- `README.md` provide curated information for understanding the codebase.
- `dev-notes.md` contains uncurated notes stored throughout development.
- Directory order indicates pipeline run order e.g. 2_* runs before 3_*.

## References

1. HypoMap3D [X]
2. [Steuernagel, L., Lam, B.Y.H., Klemm, P. et al. HypoMap—a unified single-cell gene expression atlas of the murine hypothalamus. Nat Metab 4, 1402–1419 (2022). https://doi.org/10.1038/s42255-022-00657-y]((https://www.nature.com/articles/s42255-022-00657-y))
