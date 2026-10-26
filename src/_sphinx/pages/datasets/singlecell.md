# Single Cell Objects

The HypoMap3D consists of *1,564,064* cells of which *896,959* are steady-state and *667,105* are under experimental treatment conditions. These groups are referred to as the "control" and "non-control" respectively.

<!-- ======================================================================= -->
<!--                              Full AnnData                               -->
<!-- ======================================================================= -->

## Full Single Cell AnnData

The full single cell dataset is available as an AnnData object and includes:

* Raw Gene Expression
* Hierarchical Cluster Annotations
* Latent Space Representations
* Experiment Metadata
* Gene Annotations
* UMAP Embeddings

### Structure

```plaintext
AnnData object with n_obs × n_vars = 1564064 × 56939
    X: raw counts
    obs: 'CellID', 'Dataset', 'Sample', 'total_counts', 'n_genes', 'pct_mito', 'Author', 'Source', 'Data_source_link', 'GEO', 'Publication_Link', 'Suspension', 'Technology', 'Author_Sex_harmonized', 'sex_per_sample_inferred', 'Strain', 'Diet', 'Treatment', 'Age', 'Mean_age_in_weeks', 'IsControl', 'L1_ID', 'L2_ID', 'L3_ID', 'L4_ID', 'L5_ID', 'L1_descr', 'L2_descr', 'L3_descr', 'L4_descr', 'L5_descr', 'L1_abc_like', 'L2_abc_like', 'L3_abc_like', 'L4_abc_like', 'L5_abc_like'
    var: 'gene_ids', 'IsMouseTF', 'IsHumanTFOrtholog', 'IsSourceGene', 'IsTargetGene', 'IsHvgExcludedGene', 'IsAnnotationExcludedGene', 'IsPriorityGene', 'IsTopGene', 'IsHvg'
    obsm: 'X_embed', 'X_spatial_umap', 'X_umap'
```

### Access

```{caution}
Public data access to `hm3d_sc_all.h5ad` is currently unavailable as HypoMap3D is not yet published.
```

### Using/Loading

To load the full single cell AnnData object, use the following code:

```python
import anndata as ad # or, import scanpy as sc, sc.read_h5ad()
adata = ad.read_h5ad("hm3d_sc_all.h5ad")
adata
```

### Subsetting to Control and Non-Control Cells

A boolean mask is provided in `.obs['IsControl']` to distinguish between control and non-control cells. `True` indicates a control cell, while `False` indicates a non-control cell.

For example, to subset the AnnData object to only control cells:

```python
control_adata = adata[adata.obs['IsControl']].copy()
```

And to subset to only non-control cells:

```python
non_control_adata = adata[~adata.obs['IsControl']].copy()
```

<!-- ======================================================================= -->
<!--                              Marker Genes                               -->
<!-- ======================================================================= -->

## Marker Genes Table

The marker genes for each hierarchical cluster level are provided in a CSV file, which includes the gene names and their corresponding cluster annotations.

Marker genes are computed using the Van Elteren test, a form of stratified Wilcoxon rank-sum test.

Markers were computed with respect to three cluster sets across all 5 levels independently:

1. against all other clusters of the same level (vsAll)
2. against all sister clusters of the same level (vsSisters)
3. against all cluster within the same brain region of the same level (vsRegion i.e. ZI)

### Structure

Here is the first row of the marker genes CSV file:

|clustering|celltype_id|celltype_descr|celltype_abc_like|comparison  |gene_name|specificity      |pctA              |pctB              |meanA    |meanB      |log2FC   |p  |padj|VE                |rank|
|----------|-----------|--------------|-----------------|------------|---------|-----------------|------------------|------------------|---------|-----------|---------|---|----|------------------|----|
|L1        |Gl1-1      |Gl1-1 Glut    |Gl1-1 Glut       |OneVsSisters|Foxb1    |6.236            |0.106             |0.003             |0.576    |0.005      |6.674    |0.0|0.0 |12.376            |22  |

### Access

```{caution}
Public data access to `hm3d_marker_genes_all_levels.csv` is currently unavailable as HypoMap3D is not yet published.
```

<!-- ======================================================================= -->
<!--                               Projection                                -->
<!-- ======================================================================= -->

## Projection Objects

The projection of novel data onto the HypoMap3D reference space is carried out using `scArches`. This process requires the pre-trained scVI model weights.

The projection does not require the full gene expression profiles of the reference cells. As such, we provide significantly reduced versions of the full AnnData objects, with the count matrices `.X` removed, for convenience.

### scVI Model Weights
The individual scVI model weights are required for projecting new data onto the HypoMap3D reference space using `scArches`.

The weights can be downloaded here: [XXX]().

### Reduced Full AnnData
The reduced full AnnData object is a version of the full single cell dataset with the count matrix `.X` removed. This allows for efficient storage and faster loading times while retaining all other annotations and embeddings necessary for projection.

The reduced AnnData can be downloaded here: [XXX]().

### Reduced Control AnnData
The reduced control AnnData object is a version of the control-only single cell dataset with the count matrix `.X` removed. This allows for efficient storage and faster loading times while retaining all other annotations and embeddings necessary for projection.

The reduced control AnnData can be downloaded here: [XXX]().