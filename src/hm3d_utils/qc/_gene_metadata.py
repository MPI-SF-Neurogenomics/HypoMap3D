import anndata as ad 
import pandas as pd 

def add_gene_metadata(
    adata: ad.AnnData,
    gene_metadata_df: pd.DataFrame,
):
    """
    Add gene-level metadata to an AnnData object and align it to the existing gene order.

    This function replaces ``adata.var`` with metadata from ``gene_metadata_df``,
    aligning rows by stable gene IDs and preserving the original gene order in
    ``adata``. The resulting ``adata.var`` is indexed by gene names, and
    ``adata.var_names`` is updated accordingly.

    Parameters
    ----------
    adata
        AnnData object whose ``.var`` will be updated. Must contain a
        ``gene_ids`` column in ``adata.var``.
    gene_metadata_df
        DataFrame containing gene metadata. Must include the columns
        ``"gene_ids"`` (stable gene identifiers) and ``"gene_name"``.

    Notes
    -----
    - Gene metadata is aligned using stable gene IDs, not gene names.
    - Any genes in ``adata.var`` missing from ``gene_metadata_df`` will raise
      a ``KeyError`` during alignment.
    - The index name of ``adata.var`` is reset to ``None`` for consistency.

    """

    # Add new gene metadata
    gene_metadata_mapper = gene_metadata_df.set_index("gene_ids") # use stable gene ID as index
    # Align new gene data according to previous gene order and insert into .var 
    adata.var = gene_metadata_mapper.loc[adata.var.gene_ids, :].reset_index().set_index("gene_name")
    # Update var_names
    adata.var_names = adata.var.index
    # Reset index name for comparability
    adata.var.index.name = None
