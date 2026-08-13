# Spatial Objects

Understanding the spatial localization of cells within the hypothalamus is essential both for the classification of cell types within the hypothalamus but also for understanding the interaction between them. HypoMap3D integrates both spatial and single cell transcriptomics to further explore cellular heterogeneity in the region. Raw spatial transcriptomics are combined to form an **AnnData** object and combined cluster probabilties are integrated to form a smooth voxel grids as **Zarr** objects.

<!-- ======================================================================= -->
<!--                               Raw Spatial                               -->
<!-- ======================================================================= -->

## Spatial AnnData

Single cell resolution spatial transcriptomics from 4 major datasets were concatenated and harmonized. Cells within a defined bounding box around the hypothalamus were kept leaving ~1.4 mil cells of which ~800k reside within the hypothalamus itself. The label transfer probablities of cell-types *L1-5* are included as fuzzy sparse matrices in the `.obsm`.

### Structure

```bash
AnnData object with n_obs × n_vars = 1401313 × 56939
    X: raw count matrix
    obs: 'Dataset', 'Sample', 'Author', 'Technology', 'Donor', 'Sex', 'RegionID', 'RegionAcronym', 'RegionName', 'total_counts', 'total_hvg_counts', 'total_genes', 'IsHY', 'L1_descr_max_p', 'L1_abc_like_max_p', 'L1_max_p', 'L2_descr_max_p', 'L2_abc_like_max_p', 'L2_max_p', 'L3_descr_max_p', 'L3_abc_like_max_p', 'L3_max_p', 'L4_descr_max_p', 'L4_abc_like_max_p', 'L4_max_p', 'L5_descr_max_p', 'L5_abc_like_max_p', 'L5_max_p'
    var: 'gene_stable_ids', 'IsHVG'
    uns: 'L1_ID', 'L1_abc_like', 'L1_descr', 'L2_ID', 'L2_abc_like', 'L2_descr', 'L3_ID', 'L3_abc_like', 'L3_descr', 'L4_ID', 'L4_abc_like', 'L4_descr', 'L5_ID', 'L5_abc_like', 'L5_descr', 'flags', 'pathmat'
    obsm: 'P_L1', 'P_L2', 'P_L3', 'P_L4', 'P_L5', 'X_embed', 'X_umap', 'X_umap_sc', 'coords_ccf_mm', 'coords_ccf_mm_wrap', 'coords_source', 'coords_voxels_25um', 'coords_voxels_25um_bounded'
    varm: 'IsCaptured'
```

### Access

```{caution}
Public data access is currently unavailable as dataset is not published.

```

### Using/Loading

The full anndata can be loaded using the `AnnData` python package or via the `ScanPy` wrapper

```python
import anndata as ad # via AnnData
scdata = ad.read_h5ad('sp_anndata.h5ad')

import scanpy as sc # via ScanPy
scdata = sc.read_h5ad('sp_anndata.h5ad')
```

<!-- ======================================================================= -->
<!--                              Voxel Grid(s)                              -->
<!-- ======================================================================= -->

## Voxel Grids Zarr

The smoothed spatial voxel grids are available in the Zarr file format.

### Structure

```bash
.
├── L1
│   ├── abc-like (6,) StringDType()
│   ├── descr (6,) StringDType()
│   ├── flags (6,) StringDType()
│   ├── ids (6,) StringDType()
│   ├── mass (6, 168, 144, 107) float32
│   ├── prob (6, 168, 144, 107) float32
│   └── rel_mass (6, 168, 144, 107) float32
├── L2
│   ├── abc-like (82,) StringDType()
│   ├── descr (82,) StringDType()
│   ├── flags (82,) StringDType()
│   ├── ids (82,) StringDType()
│   ├── mass (82, 168, 144, 107) float32
│   ├── prob (82, 168, 144, 107) float32
│   └── rel_mass (82, 168, 144, 107) float32
├── L3
│   ├── abc-like (274,) StringDType()
│   ├── descr (274,) StringDType()
│   ├── flags (274,) StringDType()
│   ├── ids (274,) StringDType()
│   ├── mass (274, 168, 144, 107) float32
│   ├── prob (274, 168, 144, 107) float32
│   └── rel_mass (274, 168, 144, 107) float32
├── L4
│   ├── abc-like (716,) StringDType()
│   ├── descr (716,) StringDType()
│   ├── flags (716,) StringDType()
│   ├── ids (716,) StringDType()
│   ├── mass (716, 168, 144, 107) float32
│   ├── prob (716, 168, 144, 107) float32
│   └── rel_mass (716, 168, 144, 107) float32
├── L5
│   ├── abc-like (932,) StringDType()
│   ├── descr (932,) StringDType()
│   ├── flags (932,) StringDType()
│   ├── ids (932,) StringDType()
│   ├── mass (932, 168, 144, 107) float32
│   ├── prob (932, 168, 144, 107) float32
│   └── rel_mass (932, 168, 144, 107) float32
├── bbox
│   ├── ccf_mm
│   │   ├── max (3,) float32
│   │   └── min (3,) float32
│   ├── shape (3,) float32
│   └── voxel_25_um
│       ├── max (3,) float32
│       └── min (3,) float32
├── counts (168, 144, 107) float32
├── evidence (168, 144, 107) float32
├── latent_embedding (45, 168, 144, 107) float32
└── region
    ├── acronyms (176,) StringDType()
    ├── ids (176,) int32
    ├── region_labels (168, 144, 107) uint32
    └── region_masks (176, 168, 144, 107) bool
```

### Access

```{caution}
Public data access is currently unavailable as dataset is not published. 

```

### Using/Loading

The voxel grids can be loaded using the `Zarr` python package with **version >3.1**
```python
import zarr

voxel_store = zarr.open('smoothed_voxel_grids.zarr')

# Reading slice 50 of L1 smoothed probabilities for cluster `Gb1-3 HY Gaba`
clusters = voxel_store['L1'].attrs['descr']
idx = clusters.index('Gb1-3 HY Gaba')
view = voxel_store['L1'][idx, 50, :, :]
```

```{seealso}
For more information on handling the zarr voxel grids, please refer to
the {doc}`SpatialSlicePlotting Tutorial <../tutorials/SpatialSlicePlotting>`
```
