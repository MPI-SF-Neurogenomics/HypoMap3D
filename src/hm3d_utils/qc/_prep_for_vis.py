import scanpy as sc
import anndata as ad 
from hm3d_utils.qc import compute_QC_metrics
import pandas as pd

def prepare_adata_for_visualization(adata_dict: dict, hvg_fp: None | str = None,
    highly_variable_series: None | pd.Series = None) -> ad.AnnData:
    """
    Preprocesses and prepares a concatenated AnnData object for visualization using UMAP.

    This function performs the following operations:
        1. Concatenates a dictionary of AnnData objects along a new 'Sample' index.
        2. Adds a dataset label to the `.obs` dataframe.
        3. Ensures `.var` is preserved in case it is lost during concatenation.
        4. Computes standard quality control (QC) metrics:
            - Number of genes expressed per cell
            - Total UMI counts per cell
            - Percentage of mitochondrial gene counts per cell
        5. Stores raw counts in `.layers['raw']`, or resets `.X` if raw is already present.
        6. Normalizes and log-transforms the expression data.
        7. Marks highly variable genes (HVGs) based on a precomputed list (from file or series).
        8. Performs PCA on HVGs.
        9. Constructs a KNN graph using PCA embeddings.
       10. Computes UMAP coordinates for visualization.

    Parameters
    ----------
    adata_dict : dict
        A dictionary where keys are sample identifiers and values are individual AnnData objects.
        These will be concatenated into a single object.

    hvg_fp : str or None
        File path to a text file containing a list of highly variable gene names (one per line).
        If provided, these genes will be used for downstream analysis. Must provide either this or
        `highly_variable_series`.

    highly_variable_series : pd.Series or None
        A boolean pandas Series (indexed by gene names) indicating highly variable genes.
        Must provide either this or `hvg_fp`.

    Returns
    -------
    ad.AnnData
        A single AnnData object containing:
            - Raw counts in `.layers['raw']`
            - Normalized and log-transformed data in `.X`
            - Computed QC metrics in `.obs`
            - HVG annotation in `.var['highly_variable']`
            - PCA results in `.obsm['X_pca']`, `.uns['pca']`, `.varm['PCs']`
            - KNN graph in `.obsp` and `.uns['neighbors']`
            - UMAP coordinates in `.obsm['X_umap']`
            - UMAP configuration in `.uns['umap']`

    Raises
    ------
    ValueError
        If neither `hvg_fp` nor `highly_variable_series` is provided.

    Notes
    -----
    - It is assumed that `compute_QC_metrics` is available in the global namespace.
    - This function does not perform HVG selection, only annotation based on provided input.
    """


    if (not isinstance(highly_variable_series, pd.Series)) & (not hvg_fp):
        raise ValueError("Neither hvg_pf nor highly_variable_series has been passed.")

    # Concatenate AnnData objects
    adata = sc.concat(
        adata_dict,
        #index_unique='_',
        #label='Sample',
        merge='same'
    )

    # Ensure `.var` isn't empty in case that happens during concatenation
    if adata.var.empty:
        ref_adata = next(iter(adata_dict.values()), None)
        if ref_adata is not None:
            adata.var = ref_adata.var.copy()

    # Compute quality control metrics
    compute_QC_metrics(adata)

    # Store raw counts
    if 'raw' in adata.layers:
        adata.X = adata.layers['raw'].copy()
        adata.uns.pop('log1p', None)
    else:
        adata.layers['raw'] = adata.X.copy()

    # Normalize and log-transform
    sc.pp.normalize_total(adata, target_sum=10_000)
    sc.pp.log1p(adata)

    if not hvg_fp:
        adata.var["highly_variable"] = highly_variable_series
    else:        
        # Set HVGs from file
        with open(hvg_fp) as fp:
            hvgs = [l.strip() for l in fp.readlines()]
        adata.var['highly_variable'] = adata.var_names.isin(hvgs)

    # Dimensionality reduction and neighbor graph
    sc.pp.pca(adata, n_comps=100, use_highly_variable=True)
    sc.pp.neighbors(adata, n_neighbors=15, use_rep='X_pca')
    sc.tl.umap(adata, spread=1, min_dist=0.3, random_state=0)

    return adata
