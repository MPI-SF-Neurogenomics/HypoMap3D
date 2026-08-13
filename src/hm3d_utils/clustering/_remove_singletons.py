import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
from typing import Union

def remove_singleton_clusters(
    adata: sc.AnnData,
    min_cluster_size: int = 5,
    n_iterations: int = 2,
    plot_umap: bool = True,
    umap_save_path: Union[str, Path, None] = None,
) -> sc.AnnData:
    """
    Identifies and removes singleton clusters (below a minimum size) from an AnnData object
    using Leiden clustering. Optionally generates a UMAP highlighting singleton cells.
    Assumes a UMAP is present in `adata`.

    Parameters
    ----------
    adata : sc.AnnData
        The annotated data matrix to process. This object is not modified in-place.
    min_cluster_size : int, default=5
        Clusters with fewer than this number of cells will be removed.
    n_iterations : int, default=2
        Number of iterations for the Leiden clustering algorithm.
    plot_umap : bool, default=True
        Whether to generate and display a UMAP plot highlighting singleton cells.
    umap_save_path : str or Path or None, default=None
        If provided, the path (as a string or pathlib.Path) to save the UMAP plot image.
        A ".png" extension will be automatically appended if not present.
        The parent directory will be created if it does not exist.

    Returns
    -------
    sc.AnnData
        A new AnnData object with singleton clusters removed.
    """
    # Run Leiden clustering
    sc.tl.leiden(adata, flavor='igraph', n_iterations=n_iterations)

    # Find singleton clusters
    counts = adata.obs['leiden'].value_counts()
    singleton_clusters = counts.index[counts < min_cluster_size].tolist()
    singleton_mask = adata.obs['leiden'].isin(singleton_clusters)

    # Optional: Plot UMAP showing singleton cells
    if plot_umap:
        adata.obs["IsSingleton"] = singleton_mask

        palette = {
            'True': 'red',
            'False': 'lightgrey'
        }

        fig, ax = plt.subplots(1, 1, figsize=(6, 6))
        sc.pl.umap(adata, color='IsSingleton', size=2, show=False, ax=ax, palette=palette)
        ax.set_title(f'UMAP: IsSingleton, n = {adata.obs["IsSingleton"].sum()} cells filtered')
        plt.tight_layout()
        plt.show()

    if umap_save_path is not None:
        save_path = Path(umap_save_path)
        if save_path.suffix != ".png":
            save_path = save_path.with_name(save_path.name + ".png")
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=300, bbox_inches='tight')

        adata.obs.drop(columns="IsSingleton", inplace=True)

    # Filter out singleton cells and return new object
    filtered_adata = adata[~singleton_mask].copy()
    filtered_adata.obs.drop(columns=['leiden'], inplace=True)

    return filtered_adata



