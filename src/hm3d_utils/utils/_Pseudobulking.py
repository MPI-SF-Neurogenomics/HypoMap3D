
import scanpy as sc
import pandas as pd
import numpy as np
from scipy import sparse


def create_pseudobulks(
    adata,
    sample_key="Sample",
    cell_type_key="L1",
    obs_to_keep=None
):
    """
    Create pseudobulks for each unique cell type (defined by `cell_type_key`) within each sample.

    Parameters
    ----------
    adata : AnnData
        AnnData object containing single-cell expression data.
    sample_key : str, default="Sample"
        The name of the column in `adata.obs` that contains sample or donor identifiers.
    cell_type_key : str, default="L1"
        The name of the column in `adata.obs` that defines cell type or clustering labels.
        Pseudobulks will be created per sample per cell type.
    obs_to_keep : list of str, optional
        List of metadata columns in `adata.obs` to retain in the output `.obs`.
        For each pseudobulk, the value of the first contributing cell is used.

    Returns
    -------
    AnnData
        A new AnnData object with one observation per pseudobulk (sample × cell type).
        The gene expression matrix is stored as a sparse matrix (sum of UMI counts).
        `.obs` contains:
            - sample_key (e.g., "Sample")
            - "cell_type" (the identity from `cell_type_key`)
            - "n_cells" (number of cells contributing to the pseudobulk)
            - "library_size" (sum of all gene counts)
            - any additional columns from `obs_to_keep`
            - "clustering_resolution" (value of `cell_type_key` to indicate resolution)
    """
    if obs_to_keep is None:
        obs_to_keep = []

    pseudobulk_data = []
    pseudobulk_obs = []

    # Iterate over each unique cell type
    for cell_type in adata.obs[cell_type_key].unique():
        # Subset to current cell type
        adata_subset = adata[adata.obs[cell_type_key] == cell_type].copy()
        adata_subset.obs[sample_key] = adata_subset.obs[sample_key].astype("category")

        # Aggregate within each sample
        for sample in adata_subset.obs[sample_key].cat.categories:
            adata_sample = adata_subset[adata_subset.obs[sample_key] == sample]
            if adata_sample.shape[0] == 0:
                continue

            # Sum gene expression 
            summed_expr = adata_sample.X.sum(axis=0)
            # Get library size
            library_size = summed_expr.sum()
            # Collect gene expression
            pseudobulk_data.append(sparse.csr_matrix(summed_expr))

            # Create metadata row for .obs
            obs_row = {
                sample_key: sample,
                "cell_type": cell_type,
                "n_cells": adata_sample.n_obs,
                "library_size": library_size
            }

            # Add optional metadata
            for key in obs_to_keep:
                obs_row[key] = adata_sample.obs[key].iloc[0]

            pseudobulk_obs.append(obs_row)

    # Create final AnnData object
    X_sparse = sparse.vstack(pseudobulk_data).tocsr()
    obs_df = pd.DataFrame(pseudobulk_obs)
    obs_df.index = [
        f"{cell_type_key}.{row['cell_type']}.{row[sample_key]}" for _, row in obs_df.iterrows()
    ]

    pseudobulk_adata = sc.AnnData(
        X=X_sparse,
        obs=obs_df,
        var=adata.var.copy()
    )

    # Store which resolution this pseudobulk came from
    pseudobulk_adata.obs["clustering_resolution"] = cell_type_key

    return pseudobulk_adata


def create_pseudobulks_for_multiple_levels(
    adata,
    levels=("L1", "L2", "L3", "L4", "L5"),
    sample_key="Sample",
    obs_to_keep=None
):
    """
    Create pseudobulks for each clustering level (cell type annotation column)
    and concatenate them into one AnnData object.

    Parameters
    ----------
    adata : AnnData
        Input AnnData object with single-cell data.
    levels : list or tuple
        List of column names in .obs to treat as cell type keys.
    sample_key : str
        Column in .obs denoting sample ID.
    obs_to_keep : list
        Additional metadata columns to retain in pseudobulk .obs.

    Returns
    -------
    AnnData
        Concatenated pseudobulk AnnData object across all clustering levels.
    """
    pseudobulk_list = []

    for cell_type_key in levels:
        print(f"Processing clustering level: {cell_type_key}")
        pb = create_pseudobulks(
            adata,
            sample_key=sample_key,
            cell_type_key=cell_type_key,
            obs_to_keep=obs_to_keep
        )
        pseudobulk_list.append(pb)

    pseudobulks_all = sc.concat(pseudobulk_list, axis=0)
    return pseudobulks_all

