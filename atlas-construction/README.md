# Atlas construction

This directory consists of the full HypoMap3D pipeline divided, in order, into submodules:

```bash
atlas-construction/
├── 00_setup # restoring environment, pulling containers etc
├── 01_alignment # alignment of single-cell to reference genome via CellRanger and STARsolo
├── 02_QC # Filtering of low quality and empty droplets
├── 03_SOLO # Removal of Doublets via SOLO + doublet clustering
├── 04_HY # Removal of non-hypothalamic cell-types using ABC propagation and preliminary scVI model
├── 05_clustering # Hierarchical consensus clustering and biological pruning of control cells for cell-type classification
├── 06_NC_projection # Projection of non-control into control latent space via scArches
├── 07_spatial # Projection of spatial in single-cell reference space and subsequent spatial integration
├── 08_annotation # Annotation of cell-types using regional information, ABC match, marker genes etc
├── 09_objects # Creation of final frozen data objects
└── README.md # This file
```

## Data Directory

The data directory is a folder separate from the main repository to where pipelie results are written.

For more details on the data directory and it's configuration please refere to [00_setup/README.md](./00_setup/README.md).

Each submodule constructs a corresponding results directory:
```
atlas-construction/02_QC --> $DATADIR/02_QC
```