import anndata as ad
import pandas as pd

from hm3d_utils.qc import (
    convert_to_GRCm39,
    add_gene_metadata, 
    standardize_obs_names)


def sanitize_anndata(
    adata: ad.AnnData,
    gene_metadata_df: pd.DataFrame,
    dataset_name: str,
    sample_name: str,
    map_new_genome: bool,
) -> ad.AnnData:
    """
    Sanitize and standardize an AnnData object in place.

    This function performs the following operations:
    - Ensures the count matrix is integer-typed
    - Eliminates explicit zero entries in sparse matrices
    - Optionally maps genes to a new genome build
    - Adds and aligns gene metadata in `.var`
    - Standardizes observation names and annotates `.obs`
      with CellID, Dataset, and Sample information

    Parameters
    ----------
    adata
        AnnData object to be sanitized. Modified in place.
    gene_metadata_df
        DataFrame containing gene metadata used for genome
        mapping and `.var` annotation.
    dataset_name
        Dataset identifier to store in `adata.obs["Dataset"]`.
    sample_name
        Sample identifier to store in `adata.obs["Sample"]`.
    map_new_genome
        Whether to remap genes to a new genome build.

    Returns
    -------
    The sanitized AnnData object.
    """

    # Ensure we are not modifying a view
    if adata.is_view:
        adata._init_as_actual(adata.copy())

    # Convert the count matrix to integer type
    adata.X = adata.X.astype(int)

    # Eliminate zero entries if the matrix is sparse
    if hasattr(adata.X, "eliminate_zeros"):
        adata.X.eliminate_zeros()

    # Convert to GRCm39 using gene metadata as template (in place)
    if map_new_genome:
        #print("Converting genome")
        adata = convert_to_GRCm39(adata, gene_metadata_df)
        assert adata.n_vars == len(gene_metadata_df), "Gene conversion failed: Not all .var from template are present."

    # Add / align gene metadata in .var
    add_gene_metadata(adata, gene_metadata_df)

    # Standardize obs names and add Cell/Sample/Dataset info
    standardize_obs_names(
        adata,
        dataset_name=dataset_name,
        sample_name=sample_name,
        add_info_to_obs=True,
    )

    # Clear layers to reduce size
    adata.layers.clear()    

    return adata


