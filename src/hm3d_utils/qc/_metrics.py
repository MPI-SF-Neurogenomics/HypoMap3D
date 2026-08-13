import anndata as ad 
import scanpy as sc


def compute_QC_metrics(adata: sc.AnnData) -> None:
    """
    Computes and stores quality control metrics in the AnnData object.

    Args:
        adata (sc.AnnData): The AnnData object containing single-cell RNA-seq data.

    Returns:
        None: The function modifies the AnnData object in place.
    """
    adata.var['mt'] = adata.var_names.str.lower().str.startswith('mt-')
    
    obs_qc_metrics, _ = sc.pp.calculate_qc_metrics(adata,
                                                   qc_vars='mt',
                                                   percent_top=False, 
                                                   log1p=False, 
                                                   inplace=False)
    
    adata.obs['total_counts'] = obs_qc_metrics['total_counts']
    adata.obs['n_genes'] = obs_qc_metrics['n_genes_by_counts']
    adata.obs['pct_mito'] = obs_qc_metrics['pct_counts_mt']
    
    adata.var.drop(columns='mt', inplace=True)


