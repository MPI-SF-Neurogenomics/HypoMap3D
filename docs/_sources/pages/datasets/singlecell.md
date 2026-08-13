# Single Cell Objects

The HypoMap3D consists of *1,564,064* cells of which *896,959* are steady-state and *667,105* are under experimental treatment conditions. We have divided these groups into two separate AnnData objects: control and non-control.

<!-- ======================================================================= -->
<!--                                 Control                                 -->
<!-- ======================================================================= -->

## Control AnnData

This part of the HypoMap3D dataset conists of ~897k cells classified as 'steady-state'. These cells were utilized for the definition of the cell-type taxonomy and can be utilized for exploring 'steady-state' gene dynamics and for label propagation onto new data.

### Structure

```bash
AnnData object with n_obs × n_vars = 896959 × 56939
    X: raw count matrix
    obs: 'CellID', 'Dataset', 'Sample', 'total_counts', 'n_genes', 'pct_mito', 'Author', 'Source', 'Data_source_link', 'GEO', 'Publication_Link', 'Suspension', 'Technology', 'Author_Sex_harmonized', 'sex_per_sample_inferred', 'Strain', 'Diet', 'Treatment', 'Age', 'Mean_age_in_weeks', 'IsControl', 'L1_ID', 'L2_ID', 'L3_ID', 'L4_ID', 'L5_ID', 'L1_descr', 'L2_descr', 'L3_descr', 'L4_descr', 'L5_descr', 'L1_abc_like', 'L2_abc_like', 'L3_abc_like', 'L4_abc_like', 'L5_abc_like'
    var: 'gene_ids', 'IsMouseTF', 'IsHumanTFOrtholog', 'IsSourceGene', 'IsTargetGene', 'IsHvgExcludedGene', 'IsAnnotationExcludedGene', 'IsPriorityGene', 'IsTopGene', 'IsHvg'
    obsm: 'X_embed', 'X_spatial_umap', 'X_umap'
```

### Access

```{caution}
Public data access is currently unavailable as dataset is not published. 
```

### Using/Loading

The full anndata can be loaded using the `AnnData` python package or via the `ScanPy` wrapper

```python
import anndata as ad # via AnnData
scdata = ad.read_h5ad('sc_anndata_control.h5ad')

import scanpy as sc # via ScanPy
scdata = sc.read_h5ad('sc_anndata_control.h5ad')
```

<!-- ======================================================================= -->
<!--                               Non Control                               -->
<!-- ======================================================================= -->

## Non-Control AnnData

This part of the HypoMap3D dataset conists of ~667k cells under various experimental treatments including high fat diet, influenza infection and other experiments. This provides a resource for checking potential perturbation effects from other experimental designs. Full experiment details can be found at the respective GEO and BioSampleIDs but a truncated and harmonized overview can also be found within the `anndata.obs`.

### Structure

```bash
AnnData object with n_obs × n_vars = 667105 × 56939
    X: raw count matrix
    obs: 'CellID', 'Dataset', 'Sample', 'total_counts', 'n_genes', 'pct_mito', 'Author', 'Source', 'Data_source_link', 'GEO', 'Publication_Link', 'Suspension', 'Technology', 'Author_Sex_harmonized', 'sex_per_sample_inferred', 'Strain', 'Diet', 'Treatment', 'Age', 'Mean_age_in_weeks', 'IsControl', 'L1_ID', 'L2_ID', 'L3_ID', 'L4_ID', 'L5_ID', 'L1_descr', 'L2_descr', 'L3_descr', 'L4_descr', 'L5_descr', 'L1_abc_like', 'L2_abc_like', 'L3_abc_like', 'L4_abc_like', 'L5_abc_like'
    var: 'gene_ids', 'IsMouseTF', 'IsHumanTFOrtholog', 'IsSourceGene', 'IsTargetGene', 'IsHvgExcludedGene', 'IsAnnotationExcludedGene', 'IsPriorityGene', 'IsTopGene', 'IsHvg'
    obsm: 'X_embed', 'X_spatial_umap', 'X_umap'
```

### Access

```{caution}
Public data access is currently unavailable as dataset is not published. 
```

### Using/Loading

The full anndata can be loaded using the `AnnData` python package or via the `ScanPy` wrapper

```python
import anndata as ad # via AnnData
scdata = ad.read_h5ad('sc_anndata_non_control.h5ad')

import scanpy as sc # via ScanPy
scdata = sc.read_h5ad('sc_anndata_non_control.h5ad')
```

<!-- ======================================================================= -->
<!--                                Metacell                                 -->
<!-- ======================================================================= -->

<!-- ## Metacell AnnData

The full HypoMap3D dataset is very large, as such a reduced representation was generated for quick search visualization and loading into R environments. Rather than dropping larges amount of data via subsampling, cells were grouped via KMeans to form metacells within each identified cell-type at L5. The expression of these metacells was taken as the average of the 4-80 cells they consist of.

### Structure
```bash
AnnData object with n_obs × n_vars = 40000 × 56939
    obs: 'n_cells', 'L1', 'L1_abc_like', 'L1_descr', 'L2', 'L2_abc_like', 'L2_descr', 'L3', 'L3_abc_like', 'L3_descr', 'L4', 'L4_abc_like', 'L4_descr', 'L5', 'L5_abc_like', 'L5_descr', 'pct_female', 'pct_is_nuclear'
    obsm: 'X_embed', 'X_umap'

```
```bash
[R TO BE MADE]
```

### Access
```{caution}
Public data access is currently unavailable as dataset is not published. 
```

### Using/Loading
The metacell anndata can be loaded using the `AnnData` python package or via the `ScanPy` wrapper
```python
import anndata as ad # via AnnData
scdata = ad.read_h5ad('sc_anndata_non_control.h5ad')

import scanpy as sc # via ScanPy
scdata = sc.read_h5ad('sc_anndata_non_control.h5ad')
```
The metacell anndata can be loaded in R via `Seurat`
```R
library(seurat)

...
``` -->



