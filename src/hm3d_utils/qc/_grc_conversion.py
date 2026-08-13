import scanpy as sc
import pandas as pd
import numpy as np
import multiprocessing
from functools import partial
from scipy.sparse import csr_matrix
from typing import Dict
from tqdm import tqdm



def convert_to_GRCm39(
    adata: sc.AnnData,
    template_var: pd.DataFrame,
) -> sc.AnnData:
    """
    Convert a single AnnData object to match the GRCm39 gene annotation
    provided in `template_var`.

    Uses stable gene IDs for alignment, but sets gene names as `var_names`.
    """

    # Defensive copy of var
    template_var = template_var.copy()

    # Build empty CSR matrix
    X_shape = (adata.n_obs, template_var.shape[0])
    reformatted_adata = sc.AnnData(
        X=csr_matrix(X_shape, dtype=np.float32),
        obs=adata.obs.copy(),
        var=template_var,
    )

    # Index gene IDs
    adata_gene_ids = pd.Index(adata.var["gene_ids"])
    template_gene_ids = pd.Index(template_var["gene_ids"])

    # Find overlap
    common_genes = adata_gene_ids.intersection(template_gene_ids)

    adata_idx = adata_gene_ids.get_indexer(common_genes)
    template_idx = template_gene_ids.get_indexer(common_genes)

    # Transfer expression values
    reformatted_adata[:, template_idx].X = adata[:, adata_idx].X

    # Set gene names as var_names
    reformatted_adata.var_names = reformatted_adata.var["gene_name"].astype(str)

    # Ensure uniqueness (Scanpy requirement)
    reformatted_adata.var_names_make_unique()

    return reformatted_adata


def convert_all_to_GRCm39(
    adata_dict: Dict[str, sc.AnnData],
    template_file: str,
    show_progress: bool = True
) -> Dict[str, sc.AnnData]:
    """
    Converts all AnnData objects in a dictionary to the GRCm39 genome build using a given template file.
    
    Parameters:
    -----------
    adata_dict : Dict[str, sc.AnnData]
        Dictionary of AnnData objects to convert.

    template_file : str
        Path to the reference .csv file.

    show_progress : bool
        Whether to display a progress bar.

    Returns:
    --------
    Dict[str, sc.AnnData]
        Dictionary of AnnData objects aligned to GRCm39.
    """
    # Load reference gene annotations
    template_var = pd.read_csv(template_file, index_col=0)
    template_var.feature_types= template_var.feature_types.astype("category")
    template_var.genome = template_var.genome.astype("category")

    # Partial function to fix the template_var argument
    convert_func = partial(convert_to_GRCm39, template_var=template_var)

    # Prepare input data
    keys = list(adata_dict.keys())
    adata_list = list(adata_dict.values())

    with multiprocessing.Pool() as pool:
        if show_progress:
            converted_list = list(tqdm(pool.imap(convert_func, adata_list), total=len(adata_list)))
        else:
            converted_list = pool.map(convert_func, adata_list)

    return dict(zip(keys, converted_list))
