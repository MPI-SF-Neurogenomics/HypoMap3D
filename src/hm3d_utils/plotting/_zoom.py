import warnings

import scanpy as sc
import numpy as np

from scanpy._utils import view_to_actual

from collections.abc import Sequence

def zoom_umap(adata: sc.AnnData,
              interest_mask: Sequence[bool], 
              clip_pct: float = 0.01,
              margin: float = 0.05,
              include_background=True,
              **kwargs) -> None:
    """A wrapper around `sc.pl.umap` provided a zoomed umap on a set of cells of interest.

    Parameters
    ----------
    adata : AnnData
        AnnData to plot with UMAP, must have an AnnData.obsm['X_umap'] field.
    interest_mask : np.array[bool]
        A boolean mask of cells of interest
    clip_pct : float, optional
        Percentile of cells to ignore when defining bounding box, by default 0.01
    margin : float, optional
        Added padding around bounding as a percentage of zoom span, by default 0.05
    include_background: bool, optional
        Whether to plot greyed out points not part of interest mask but within bounding box, by default True
    kwargs : any, optional
        Keyword arguments passed to `sc.pl.umap`.
    """

    X_umap = adata.obsm['X_umap']

    # Define zoom bounding box --------------------------------------------------- #
    bbox_min = np.quantile(X_umap[interest_mask], clip_pct, axis=0)
    bbox_max = np.quantile(X_umap[interest_mask], 1-clip_pct, axis=0)

    bbox_span = bbox_max - bbox_min

    bbox_min = bbox_min - margin * bbox_span # a border around min max cells
    bbox_max = bbox_max + margin * bbox_span

    # Subset AnnData to inbounds cells ------------------------------------------- #
    inbounds_mask = (X_umap >= bbox_min).all(axis=1) & (X_umap <= bbox_max).all(axis=1)

    # subset = adata[inbounds_mask] # could call copy...
    # with warnings.catch_warnings():
    #     warnings.simplefilter("ignore")
    #     view_to_actual(subset)
    inbound_interest_mask = interest_mask & inbounds_mask
    inbound_background_mask = ~interest_mask & inbounds_mask

    # Plot via sc.pl.umap -------------------------------------------------------- #

    if include_background:
        bg_kwargs = kwargs.copy()
        bg_kwargs.pop('color')
        bg_kwargs.pop('add_outline')
        # bg_kwargs['color'] = "lightgrey"
        # bg_kwargs['add_outline'] = False
        sc.pl.umap(adata[inbound_background_mask].copy(), **bg_kwargs)

    sc.pl.umap(adata[inbound_interest_mask].copy(), **kwargs)