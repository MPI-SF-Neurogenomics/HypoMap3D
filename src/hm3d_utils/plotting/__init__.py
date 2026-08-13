from ._colors import cmap_marker
from ._contour import plot_clusters_collapsed, plot_single_cluster
from ._embedding_labels import embedding_adjustedlabels
from ._embedding_zoom import zoom_embedding, zoom_embedding_circular

from ._qc_plots import (
    plot_qc_umaps,
    plot_qc_filtered_vs_clean,
    plot_qc_filtered_vs_manual_vs_clean,
    plot_qc_histogram,
    plot_qc_thr_hists
)

# from ._label_transfer import plot_label_propagation_summary

from ._zoom import zoom_umap

__all__ = [
    cmap_marker,
    plot_clusters_collapsed,
    plot_single_cluster,
    embedding_adjustedlabels,
    zoom_embedding, zoom_embedding_circular,
    plot_qc_umaps,
    plot_qc_filtered_vs_clean,
    plot_qc_filtered_vs_manual_vs_clean,
    plot_qc_histogram,
    plot_qc_thr_hists,
    # plot_label_propagation_summary,
    zoom_umap
]