import numpy as np
import anndata as ad

def clean_and_export_adata(
    adata: ad.AnnData,
    min_counts_threshold: int,
    min_genes_threshold: int,
    max_pct_mito_threshold: float,
    manually_removed_obs: list
) -> ad.AnnData:
    """
    Clean an AnnData object by filtering out failed cells, trimming unnecessary attributes,
    and preparing for downstream analysis or export.

    Parameters:
    - adata: AnnData object to be cleaned
    - min_counts_threshold: Minimum UMI counts per cell
    - min_genes_threshold: Minimum number of genes per cell
    - max_pct_mito_threshold: Maximum mitochondrial content allowed
    - manually_removed_obs: List of manually flagged cells to record in metadata

    Returns:
    - export_adata: Cleaned AnnData object
    """
    # Mark all cells to be removed (failed QC or manual)
    adata.obs["failed_qc"] = adata.obs["failed_qc"] | adata.obs["manually_removed"]

    # Subset the data to only good-quality cells
    export_adata = adata[~adata.obs["failed_qc"]].copy()
    export_adata.obs.drop(columns="failed_qc", inplace=True)

    # remove leiden clusterings
    export_adata.obs = export_adata.obs.loc[:, [c for c in list(export_adata.obs.columns) if "leiden" not in c]]

    # Remove unused fields in .uns
    export_adata.uns = {
        k: v for k, v in export_adata.uns.items()
        if ("color" not in k and k != "log1p" and "pca" not in k and "neighbors" not in k and "leiden" not in k)
    }

    # Remove unused fields in .obsm
    export_adata.obsm = {
        k: v for k, v in export_adata.obsm.items()
        if ("pca" not in k)
    }

    # get rid of fake highly variable genes
    export_adata.var = export_adata.var.loc[:, [c for c in export_adata.var.columns if "variable" not in c]]

    # delete unused fields
    del export_adata.obsp, export_adata.varm

    # Set .X to raw counts (copied from .layers)
    export_adata.X = export_adata.layers["raw"].copy()

    # Ensure .X is integer-typed
    if not np.issubdtype(export_adata.X.data.dtype, np.integer):
        export_adata.X = export_adata.X.astype(np.int32)

    # Clear layers to reduce size
    export_adata.layers.clear()

    # Add thresholds and manual removal record to .uns
    export_adata.uns["min_counts_threshold"] = min_counts_threshold
    export_adata.uns["min_genes_threshold"] = min_genes_threshold
    export_adata.uns["max_pct_mito_threshold"] = max_pct_mito_threshold
    export_adata.uns["manually_removed_obs"] = manually_removed_obs

    return export_adata
