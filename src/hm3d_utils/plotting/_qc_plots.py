import os
import matplotlib.pyplot as plt
import scanpy as sc
import matplotlib.ticker as ticker
import numpy as np

def plot_qc_umaps(adata, outdir):
    """
    Plots UMAPs of key QC metrics from an AnnData object and saves the figure.

    Parameters:
    - adata: AnnData object
    - outdir: Path to the output directory (figure will be saved as outdir + '/umap_qc_metrics.png')
    """
    # Create output path if it doesn't exist
    os.makedirs(outdir, exist_ok=True)
    outfile = os.path.join(outdir, "umap_qc_metrics.png")

    # Create subplots
    fig, axes = plt.subplots(1, 4, figsize=(20, 5))

    # Plot each QC metric
    sc.pl.umap(adata, color='total_counts', ax=axes[0], show=False, title="total_counts", vmax=4000, cmap='viridis')
    sc.pl.umap(adata, color='n_genes', ax=axes[1], show=False, title="n_genes", vmax=3000, cmap='inferno')
    sc.pl.umap(adata, color='pct_mito', ax=axes[2], show=False, title="pct_mito", vmax=20, cmap='YlOrRd')
    sc.pl.umap(adata, color='Sample', ax=axes[3], show=False, title="SampleID", vmax=20, palette='Dark2', legend_loc=None)

    # Adjust layout and save/show figure
    plt.tight_layout()
    plt.savefig(outfile, dpi=300, bbox_inches='tight')
    plt.show()

    print(f"UMAP QC metrics saved to: {outfile}")


def plot_qc_filtered_vs_clean(adata, outdir):
    """
    Plots UMAPs showing QC filtering status and clean data by sample.
    
    Parameters:
    - adata: AnnData object (must have 'failed_qc' in .obs)
    - outdir: Path to the output directory (figure saved as outdir + '/umap_filtered_vs_clean.png')
    """
    # Create output path if it doesn't exist
    os.makedirs(outdir, exist_ok=True)
    outfile = os.path.join(outdir, "umap_filtered_vs_clean.png")

    # Define custom palette for QC filtering
    selected_palette = {
        'True': 'red',
        'False': 'lightgrey'
    }

    # Create subplots
    fig, axes = plt.subplots(1, 2, figsize=(14, 7))

    # Plot filtered view
    sc.pl.umap(
        adata,
        color='failed_qc',
        ax=axes[0],
        show=False,
        title="filtered",
        palette=selected_palette,
        legend_loc=None
    )

    # Plot clean cells only
    sc.pl.umap(
        adata[~adata.obs['failed_qc']],
        color='Sample',
        ax=axes[1],
        show=False,
        title="clean",
        palette='Dark2'
    )

    # Save and show
    plt.tight_layout()
    plt.savefig(outfile, dpi=300, bbox_inches='tight')
    plt.show()

    print(f"Filtered vs Clean UMAP saved to: {outfile}")


def plot_qc_filtered_vs_manual_vs_clean(adata, outdir, 
    clusters_to_remove:list | None = None, 
    removal_clustering_res: int | float | None = None ):
    """
    Generate a 3-panel UMAP plot comparing cells removed by default QC filtering,
    manual cluster-based removal, and the final cleaned dataset after both filters.

    Parameters:
    ----------
    adata : AnnData
        Annotated data matrix. Must include boolean columns 'failed_qc' and 
        'manually_removed' in `adata.obs` indicating cells flagged by each method.
    
    outdir : str
        Path to the directory where the UMAP figure will be saved. The plot will be
        saved as 'umap_filtered_vs_manual_vs_clean.png' inside this directory.
    
    clusters_to_remove : list, optional
        List of cluster identifiers that were manually selected for removal.
        If provided, this list will be shown in the title of the second UMAP panel.

    removal_clustering_res : int or float, optional
        The resolution value used for clustering when selecting clusters to remove.
        Only used for title annotation if `clusters_to_remove` is provided.

    Returns:
    -------
    None
        The function saves the figure and displays it, but does not return any value.
    """
    # Create output path
    os.makedirs(outdir, exist_ok=True)
    outfile = os.path.join(outdir, "umap_filtered_vs_manual_vs_clean.png")

    # Define the palette
    selected_palette = {
        'True': 'red',
        'False': 'lightgrey'
    }

    # Create the figure
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # Plot 1: Original QC filter
    sc.pl.umap(
        adata,
        color='failed_qc',
        ax=axes[0],
        show=False,
        title="filtered",
        palette=selected_palette,
        legend_loc=None
    )

    # Plot 2: Manual removals
    if clusters_to_remove != None:
        title2 = f"removed manually: \nclusters {clusters_to_remove} from resolution = {removal_clustering_res}"
    else:
        title2 = "removed manually"
    sc.pl.umap(
        adata,
        color='manually_removed',
        ax=axes[1],
        show=False,
        title=title2,
        palette=selected_palette,
        legend_loc=None
    )

    # Plot 3: Clean after combining both filters
    combined_filter = adata.obs['failed_qc'] | adata.obs['manually_removed']
    sc.pl.umap(
        adata[~combined_filter],
        ax=axes[2],
        show=False,
        title="clean"
    )

    # Final layout & save
    plt.tight_layout()
    plt.savefig(outfile, dpi=300, bbox_inches='tight')
    plt.show()

    print(f"UMAP filtered vs manual vs clean saved to: {outfile}")

def plot_qc_thr_hists(
    df, 
    bins=100, 
    min_counts_threshold=None, 
    min_genes_threshold=None, 
    max_pct_mito_threshold=None, 
    save_dir=None,
    dataset_id=None  
):
    """
    Plots histograms for 'total_counts', 'n_genes', and 'pct_mito' with log10 scaling.
    
    Parameters:
        df (pd.DataFrame): Must contain 'total_counts', 'n_genes', 'pct_mito'.
        bins (int): Number of bins in the histograms.
        min_counts_threshold (float or None): Optional threshold for 'total_counts'.
        min_genes_threshold (float or None): Optional threshold for 'n_genes'.
        max_pct_mito_threshold (float or None): Optional threshold for 'pct_mito'.
        save_dir (str or None): Directory path to save the figure. If None, the figure is only shown.
        dataset_id (str or None): Identifier for the dataset, shown as the figure title.
    """
    features = [
        ('total_counts', min_counts_threshold, True),
        ('n_genes', min_genes_threshold, True),
        ('pct_mito', max_pct_mito_threshold, False)
    ]
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    for ax, (feature, threshold, log_x) in zip(axes, features):
        plot_qc_histogram(ax, df[feature], feature, threshold, bins, log_x)

    # Add suptitle
    if dataset_id:
        fig.suptitle(f"QC Threshold Histograms | Dataset: {dataset_id}", fontsize=16, y=1.05)
    
    plt.tight_layout()

    # Check if save_dir exists, if not, create it
    if save_dir:
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
            print(f"Directory {save_dir} created.")
        
        # Save figure
        save_path = os.path.join(save_dir, "qc_thr.png")
        plt.savefig(save_path, bbox_inches='tight')
        print(f"Figure saved to {save_path}")
    
    # Always show the plot
    plt.show()



def plot_qc_histogram(ax, values, feature, threshold, bins, log_x):
    """
    Helper function to plot a histogram on a given axis with appropriate scaling.
    
    Parameters:
        ax (matplotlib.axes.Axes): Axis to plot on.
        values (pd.Series or np.ndarray): Data values to plot.
        feature (str): Feature name for axis labels and title.
        threshold (float or None): Optional threshold line.
        bins (int): Number of bins in the histogram.
        log_x (bool): Whether to use a log scale for the x-axis.
    """
    values = values[values > 0]  # Avoid log of zero or negative values

    # Use log-spaced bins if log_x is True
    if log_x:
        min_val = values.min()
        max_val = values.max()
        bin_edges = np.logspace(np.log10(min_val), np.log10(max_val), bins + 1)
    else:
        bin_edges = np.linspace(values.min(), values.max(), bins + 1)

    counts, _ = np.histogram(values, bins=bin_edges)

    # Plot histogram bars with uniform bin widths in the current scale
    ax.bar(bin_edges[:-1], counts, width=np.diff(bin_edges), align='edge', edgecolor='black')

    # Apply axis scales
    if log_x:
        ax.set_xscale('log')
        ax.xaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=15))
    
    ax.set_yscale('log')
    ax.yaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=15))

    # Labeling
    ax.set_ylabel("Count")
    ax.set_xlabel(f"{feature}")
    ax.set_title(f"{feature} | thr: {threshold}" if threshold is not None else feature)

    # Threshold line
    if threshold is not None and threshold > 0:
        ax.axvline(threshold, color='red', linestyle='-', linewidth=2)
