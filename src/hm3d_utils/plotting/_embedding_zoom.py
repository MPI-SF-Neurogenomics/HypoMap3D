import scanpy as sc
import numpy as np

from collections.abc import Sequence


def zoom_embedding_circular(
    adata: sc.AnnData,
    basis: str,
    interest_mask: Sequence[bool],
    clip_pct: float = 0.01,
    margin: float = 0.05,
    include_background: bool = True,
    **kwargs
) -> None:

    X_umap = adata.obsm['X_umap']
    X_interest = X_umap[interest_mask]

    # ------------------------------------------------------------
    # Define circular region
    # ------------------------------------------------------------

    # Center: centroid of interest cells
    center = np.median(X_interest, axis=0)

    # Distances from centroid
    dists = np.linalg.norm(X_interest - center, axis=1)

    # Radius from high percentile of distances
    radius = np.quantile(dists, 1 - clip_pct)

    # Add margin
    radius *= (1 + margin)

    # ------------------------------------------------------------
    # Identify in-circle cells
    # ------------------------------------------------------------

    # Distances for all cells
    all_dists = np.linalg.norm(X_umap - center, axis=1)

    in_circle = all_dists <= radius
    inbound_interest_mask = interest_mask & in_circle
    inbound_background_mask = (~interest_mask) & in_circle

    # ------------------------------------------------------------
    # Plot
    # ------------------------------------------------------------

    
    # We need to pass ax explicitly to Scanpy

    if include_background and np.any(inbound_background_mask):
        bg_kwargs = kwargs.copy()
        bg_kwargs.pop("color", None)
        bg_kwargs.pop("add_outline", None)
        bg_kwargs.pop('show', None)
        bg_kwargs.setdefault('size', 120000 / adata.shape[0])

        sc.pl.embedding(adata[inbound_background_mask].copy(), basis=basis, show=False, **bg_kwargs)

    sc.pl.embedding(adata[inbound_interest_mask].copy(), basis=basis, **kwargs)


def zoom_embedding(adata: sc.AnnData,
                   basis: str,
                   interest_mask: Sequence[bool] | None = None, 
                   clip_pct: float = 0.01,
                   margin: float = 0.05,
                   include_background: bool = True,
                   **kwargs) -> None:

    X_embedding = adata.obsm[basis]

    if interest_mask is None:
        interest_mask = np.ones(adata.shape[0]).astype(bool)

    # Define zoom bounding box --------------------------------------------------- #
    bbox_min = np.quantile(X_embedding[interest_mask], clip_pct, axis=0)
    bbox_max = np.quantile(X_embedding[interest_mask], 1-clip_pct, axis=0)

    bbox_span = bbox_max - bbox_min

    bbox_min = bbox_min - margin * bbox_span # a border around min max cells
    bbox_max = bbox_max + margin * bbox_span

    # Subset AnnData to inbounds cells ------------------------------------------- #
    inbounds_mask = (X_embedding >= bbox_min).all(axis=1) & (X_embedding <= bbox_max).all(axis=1)

    # subset = adata[inbounds_mask] # could call copy...
    # with warnings.catch_warnings():
    #     warnings.simplefilter("ignore")
    #     view_to_actual(subset)
    inbound_interest_mask = interest_mask & inbounds_mask
    inbound_background_mask = ~interest_mask & inbounds_mask

    # Plot via sc.pl.umap -------------------------------------------------------- #

    if include_background and (np.sum(inbound_background_mask) > 0):
        bg_kwargs = kwargs.copy()
        bg_kwargs.pop('color', None)
        bg_kwargs.pop('add_outline', None)
        bg_kwargs.pop('show', None)
        bg_kwargs.setdefault('size', 120000 / adata.shape[0])

        sc.pl.embedding(adata[inbound_background_mask].copy(), basis=basis, show=False, **bg_kwargs)

    sc.pl.embedding(adata[inbound_interest_mask].copy(), basis=basis, **kwargs)