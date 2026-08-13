from typing import Optional
import pandas as pd
import scanpy as sc
from anndata import AnnData

from hm3d_utils.clustering import propagate_labels_within_dataset

def get_neuron_biased_hvgs(
    adata: AnnData,
    n_unspec: int,
    n_neuron: int,
    outfile: str,
) -> pd.DataFrame:
    """
    Identify neuron-biased highly variable genes (HVGs) from a single-cell dataset,
    combining unspecific and neuron-specific HVGs with exclusion of predetermined genes to be excluded in .var.

    To distinguish between cell types we propagate ABC class labels from the Yao dataset to all other cells, 
    and use those labels to subset the data.

    Parameters:
        adata (AnnData): Annotated data matrix.
        n_unspec (int): Number of unspecific HVGs to select.
        n_neuron (int): Number of neuron-specific HVGs to select.
        outfile (str): Path to save the final gene list.

    Returns:
        pd.DataFrame: DataFrame containing selected HVGs.
    """
    propagation_label = 'class_id_label'
    propagation_key = 'ABC_class'

    # Label propagation step 
    propagate_labels_within_dataset(
        adata, propagation_label=propagation_label, propagation_key=propagation_key
    )

    # Define non-neuronal (glial) labels
    glia_labels = ['30 Astro-Epen', '31 OPC-Oligo', '33 Vascular', '34 Immune']

    # Annotate neuronal cells
    adata.obs["IsNeuronal"] = ~(
        adata.obs[propagation_label].isin(glia_labels) |
        adata.obs[f"{propagation_key}_propagated"].isin(glia_labels)
    )

    # Filter out samples with fewer than 2 cells
    sample_counts = adata.obs['Sample'].value_counts()
    valid_samples = sample_counts[sample_counts >= 2].index
    adata = adata[adata.obs['Sample'].isin(valid_samples)].copy()

    # Compute highly variable genes (HVGs) across all cell types
    hvg_unspec_df = sc.pp.highly_variable_genes(
        adata,
        inplace=False,
        n_top_genes=adata.shape[1],
        batch_key='Sample',
        flavor='seurat_v3',
        span=1
    ).sort_values('highly_variable_rank')

    # Subset to neuronal cells only
    adata_neuron = adata[adata.obs["IsNeuronal"]].copy()

    # Filter again for neuronal cells by sample count
    neuron_sample_counts = adata_neuron.obs['Sample'].value_counts()
    valid_neuron_samples = neuron_sample_counts[neuron_sample_counts >= 2].index
    adata_neuron = adata_neuron[adata_neuron.obs['Sample'].isin(valid_neuron_samples)].copy()

    # Compute HVGs for neuronal subset
    hvg_neuron_df = sc.pp.highly_variable_genes(
        adata_neuron,
        inplace=False,
        n_top_genes=adata_neuron.shape[1],
        batch_key='Sample',
        flavor='seurat_v3',
        span=1
    ).sort_values('highly_variable_rank')

    # Exclude specified HVGs
    excluded_genes = list(adata.var[adata.var["IsHvgExcludedGene"]].index)

    # Remove excluded genes from unspecific HVGs
    hvg_unspec_df.loc[hvg_unspec_df.index.isin(excluded_genes), 'highly_variable_rank'] = None
    hvg_unspec_df = hvg_unspec_df.sort_values('highly_variable_rank')

    # Remove excluded genes from neuron HVGs
    hvg_neuron_df.loc[hvg_neuron_df.index.isin(excluded_genes), 'highly_variable_rank'] = None
    hvg_neuron_df = hvg_neuron_df.sort_values('highly_variable_rank')

    # Select top HVGs
    top_unspec = hvg_unspec_df.head(n_unspec)
    top_neuron = hvg_neuron_df.loc[~hvg_neuron_df.index.isin(top_unspec.index)].head(n_neuron)

    # Concatenate and trim to expected size
    final_df = pd.concat([top_unspec, top_neuron]).head(n_unspec + n_neuron)

    # Save gene names to file
    final_df['gene_name'].to_csv(outfile, index=False, header=False)

    return final_df
