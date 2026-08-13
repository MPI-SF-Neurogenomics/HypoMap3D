import pandas as pd
import anndata as ad
from typing import Literal
from ._markers import find_markers
from ..utils._ClusterTreeNode import ClusterTree, combine_cluster_nodes

def get_final_OneVsSister_markers_(adata: ad.AnnData, lookup: ClusterTree,
                                  cell_type: Literal["neuronal", "nonneuronal"], output_file: str, 
                                  max_clustering_level: int = 5, 
                                  strata_key: str = 'Sample',
                                  max_stratum_size: int = 1000,
                                  min_stratum_size: int = 10,
                                  min_pct_expression: float = 0.1,
                                  min_specificity: float = 0.1,
                                  padj_threshold: float = 0.05,
                                  log2fc_threshold: float = 0.0,
                                  write_header: bool = True):
    """
    Identify and save marker genes by comparing each cluster to its sister clusters 
    (One-vs-Sisters) in a hierarchical clustering tree.

    This function iterates over clusters in the provided cluster tree, performs differential
    expression marker analysis comparing each cluster to its sister clusters using specified
    thresholds and parameters, and appends the significant marker results to a CSV file.

    Parameters
    ----------
    adata : ad.AnnData
        Annotated data matrix containing gene expression data.
    lookup : ClusterTree
        Cluster tree object containing clusters with hierarchical relationships.
    cell_type : {"neuronal", "nonneuronal"}
        Subset of cells to analyze: 'neuronal' or 'nonneuronal'.
    output_file : str
        Path to the CSV file where marker gene results will be saved (appended).
    max_clustering_level : int, optional
        Maximum hierarchical clustering level to process clusters from (default is 5).
    strata_key : str, optional
        Key in `adata.obs` used to stratify samples during marker analysis (default is 'Sample').
    max_stratum_size : int, optional
        Maximum number of cells per stratum during marker detection (default is 1000).
    min_stratum_size : int, optional
        Minimum number of cells per stratum to include in marker analysis (default is 10).
    min_pct_expression : float, optional
        Minimum fraction of cells expressing a gene to consider it (default is 0.1).
    min_specificity : float, optional
        Minimum specificity score required for markers (default is 0.1).
    padj_threshold : float, optional
        Adjusted p-value cutoff to determine significance (default is 0.05).
    log2fc_threshold : float, optional
        Log2 fold-change threshold for markers to be considered upregulated (default is 0.0).
    write_header : bool, optional
        Whether to write the CSV header (default is True). Header is written only once.

    Returns
    -------
    None
        Results are written to the specified CSV file. No value is returned.

    Notes
    -----
    - Clusters above `max_clustering_level` are skipped.
    - Clusters without sister clusters are skipped.
    - If multiple sisters exist, they are combined into a single group for comparison.
    - Only significant markers passing thresholds for adjusted p-value and fold change
      are saved.
    """

    # restrict anndata object to the cell type in question such that the .indexer returned later match the index of the adata
    nn_mask = adata.obs['L1'].str.contains("NN")
    if cell_type == "neuronal":
        cell_type_adata = adata[~nn_mask]
    else:
        cell_type_adata = adata[nn_mask]

    for name, cluster in lookup.items():

            print(f"Processing cluster {name}")

            clustering_level = cluster.level
            # Dont process the cluster if not warranted.
            if clustering_level > max_clustering_level:
                print(f"Cluster {name} is from tree level {clustering_level} which is > {max_clustering_level}. Skipping.")
                continue

            else:
                sisters = cluster.get_sisters()

                if not sisters:
                    print(f"Skipping cluster {name} — no sisters found")
                    continue

                # Combine sisters if more than one
                if len(sisters) > 1:
                    sisters = combine_cluster_nodes(sisters, 'sisters')
                else:
                    sisters = sisters[0]

                ranking = find_markers(
                    cell_type_adata[cluster.indexer],
                    cell_type_adata[sisters.indexer],
                    strata_key=strata_key,
                    max_stratum_size=max_stratum_size,
                    min_stratum_size=min_stratum_size,
                    min_pct_expression=min_pct_expression,
                    min_specificity=min_specificity
                )

                # Validate ranking result
                if (
                    isinstance(ranking, pd.DataFrame) and
                    not ranking.empty and
                    {'padj', 'log2FC', 'specificity'}.issubset(ranking.columns)
                ):
                    # Filter significant, upregulated markers
                    ranking = ranking[
                        (ranking['padj'] < padj_threshold) &
                        (ranking['log2FC'] > log2fc_threshold)
                    ].sort_values('specificity', ascending=False)

                    if not ranking.empty:
                        ranking['clustering'] = f"L{clustering_level}"
                        ranking['cluster'] = cluster.name
                        ranking["comparison"] = "OneVsSisters"

                        # Append to output CSV
                        ranking.to_csv(
                            output_file,
                            mode='a',
                            header=write_header,
                            index_label='gene_name'
                        )

                        write_header = False  # Only write header on first write



def get_final_OneVsAll_markers_(
    adata: ad.AnnData,
    cell_type: Literal["neuronal", "nonneuronal"],
    output_file: str,
    max_clustering_level: int = 5,
    strata_key: str = 'Sample',
    max_stratum_size: int = 1000,
    min_stratum_size: int = 10,
    min_pct_expression: float = 0.1,
    min_specificity: float = 0.1,
    padj_threshold: float = 0.05,
    log2fc_threshold: float = 0.0,
    write_header: bool = True
):
    """
    Identify and save marker genes for each cluster against all other cells 
    within specified neuronal or non-neuronal subsets.

    This function subsets the AnnData object by cell type, iterates through hierarchical
    clustering levels, performs differential marker analysis for each cluster versus all 
    other cells, and writes significant markers to a CSV file.

    Parameters
    ----------
    adata : ad.AnnData
        Annotated data matrix containing gene expression and metadata.
    cell_type : {"neuronal", "nonneuronal"}
        Subset of cells to analyze: 'neuronal' or 'nonneuronal'.
    output_file : str
        File path to append the marker results CSV.
    max_clustering_level : int, optional
        Maximum clustering hierarchy level to process (default is 5).
    strata_key : str, optional
        Observation key in `adata.obs` used for stratifying samples (default 'Sample').
    max_stratum_size : int, optional
        Maximum number of cells per stratum during marker testing (default 1000).
    min_stratum_size : int, optional
        Minimum number of cells per stratum to include (default 10).
    min_pct_expression : float, optional
        Minimum fraction of cells expressing a gene to consider it (default 0.1).
    min_specificity : float, optional
        Minimum specificity score to filter markers (default 0.1).
    padj_threshold : float, optional
        Adjusted p-value cutoff for significance (default 0.05).
    log2fc_threshold : float, optional
        Minimum log2 fold-change threshold for marker genes (default 0.0).
    write_header : bool, optional
        Whether to write the CSV header (default True). Header is only written once.

    Returns
    -------
    None
        Writes marker gene results to the specified CSV file.
    
    Raises
    ------
    ValueError
        If 'L1' column is missing from `adata.obs`, necessary for cell type subsetting.

    Notes
    -----
    - Cells are subsetted based on presence of "NN" in `adata.obs['L1']` to distinguish
      neuronal and non-neuronal cells.
    - Marker detection is run for each cluster at each clustering level.
    - Only markers passing adjusted p-value and log2 fold-change thresholds are saved.
    """


    # subset the data to the neuronal/non-neuronal cells
    if 'L1' not in adata.obs.columns:
        raise ValueError("Expected 'L1' column in adata.obs for cell type subsetting")

    nn_mask = adata.obs['L1'].str.contains("NN")
    if cell_type == "neuronal":
        cell_type_adata = adata[~nn_mask]
    else:
        cell_type_adata = adata[nn_mask]

    clustering_df = cell_type_adata.obs.loc[:, [f"L{i}" for i in range(1,max_clustering_level + 1)]]

    # Loop over all clustering columns
    for clustering in clustering_df.columns:

        print(f"Processing clustering: {clustering}")
        # Get unique labels
        cluster_labels = clustering_df[clustering].unique()

        for cluster in cluster_labels:
            print(f"Processing cluster: {cluster}")

            mask = cell_type_adata.obs[clustering] == cluster
            cluster_data = cell_type_adata[mask]
            other_data = cell_type_adata[~mask]

            # Run marker test
            ranking = find_markers(
                cluster_data,
                other_data,
                strata_key=strata_key,
                max_stratum_size=max_stratum_size,
                min_stratum_size=min_stratum_size,
                min_pct_expression=min_pct_expression,
                min_specificity=min_specificity
            )

            # Check ranking is a non-empty DataFrame with required columns
            if (
                isinstance(ranking, pd.DataFrame) and 
                not ranking.empty and 
                {'padj', 'log2FC', 'specificity'}.issubset(ranking.columns)
            ):
                # Filter significant, upregulated genes
                ranking = ranking[
                    (ranking['padj'] < padj_threshold) &
                    (ranking['log2FC'] > log2fc_threshold)
                ].sort_values('specificity', ascending=False)

                if not ranking.empty:
                    # Add cluster info
                    ranking['clustering'] = clustering
                    ranking['cluster'] = cluster
                    ranking["comparison"] = "OneVsAll"

                    # Write to CSV
                    ranking.to_csv(
                        output_file,
                        mode='a',
                        header=write_header,
                        index_label='gene_name'
                    )

                    write_header = False  # Only write header the first time



def get_final_markers(
    adata: ad.AnnData, 
    n_lookup: ClusterTree, 
    nn_lookup: ClusterTree, 
    output_file: str
):
    """
    Run marker detection for both neuronal and non-neuronal clusters and append results to a CSV file.

    This function runs two types of marker analyses One-vs-Sister clusters and One-vs-All cells—
    separately for neuronal and non-neuronal subsets of the data.

    Parameters
    ----------
    adata : ad.AnnData
        Annotated data matrix with expression and metadata.
    n_lookup : ClusterTree
        Cluster lookup/tree for neuronal cells.
    nn_lookup : ClusterTree
        Cluster lookup/tree for non-neuronal cells.
    output_file : str
        Path to the output CSV file where marker results are appended.

    Returns
    -------
    None
        Writes marker gene results to the specified CSV file.
    """


    print("Processing neuronal clusters")

    get_final_OneVsSister_markers_(
        adata = adata, 
        lookup = n_lookup, 
        cell_type = "neuronal",
        output_file = output_file, 
        max_clustering_level = 5, 
        write_header = True # only write the header when the file is created
        )

    get_final_OneVsAll_markers_(
        adata = adata, 
        cell_type = "neuronal",
        output_file = output_file,
        max_clustering_level = 5,
        write_header = False # only write the header when the file is created
    )

    print("Processing non-neuronal clusters")
    get_final_OneVsSister_markers_(
        adata = adata, 
        lookup = nn_lookup, 
        cell_type = "nonneuronal",
        output_file = output_file, 
        max_clustering_level = 4, 
        write_header = False # only write the header when the file is created
        )

    get_final_OneVsAll_markers_(
        adata = adata, 
        cell_type = "nonneuronal",
        output_file = output_file,
        max_clustering_level = 4,
        write_header = False # only write the header when the file is created
    )
    
