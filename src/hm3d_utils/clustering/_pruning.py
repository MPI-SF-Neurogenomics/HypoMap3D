import numpy as np

from dataclasses import dataclass
from sklearn.neighbors import NearestNeighbors
from scipy.sparse import csr_matrix

from ..statistics import van_elteren_test

# ============================================================================ #
#                               Exposed Class(es)                              #
# ============================================================================ #

@dataclass
class PruningCluster:
    """ [ChatGPT]
    Lightweight container representing a cluster during iterative pruning.

    Each instance tracks both geometric (embedding centroid) and expression
    statistics used to determine whether two clusters should be merged.

    The object is intentionally immutable in spirit: merging always produces
    a new PruningCluster rather than modifying existing ones.

    Attributes
    ----------
    name : str
        Unique identifier for the cluster (used to track tested merge pairs).
    indices : np.ndarray
        Indices of observations (cells) belonging to the cluster.
    centroid : np.ndarray
        Mean embedding position of the cluster (e.g. UMAP / PCA coordinates).
    pct_gexpr : np.ndarray
        Fraction of cells expressing each gene (>0).
    mean_gexpr : np.ndarray
        Mean gene expression across cells (not log-transformed).

    Notes
    -----
    These statistics are maintained explicitly to avoid repeated sparse
    slicing during merge testing, which is the computational bottleneck.
    """
    name: str
    indices: np.ndarray
    centroid: np.ndarray
    pct_gexpr: np.ndarray
    mean_gexpr: np.ndarray
    
    @property
    def size(self):
        return len(self.indices)


# ============================================================================ #
#                              Exposed Function(s)                             #
# ============================================================================ #

def prune_clusters(
    membership: np.ndarray,
    gexpr: csr_matrix,
    embedding: np.ndarray,
    strata: np.ndarray,
    nn_checks: int,
    min_cluster_size: int,
    max_stratum_size: int,
    min_stratum_size: int,
    min_strata: int,
    min_pct_expression: float,
    min_specificity: float,
    significance_threshold: float,
    min_marker_genes: int,
    symmetric: bool = True,
    random_state: int = 0,
) -> np.ndarray:
    """ [ChatGPT]
    Iteratively merges transcriptionally indistinguishable clusters using
    local neighborhood search and stratified differential expression testing.

    The algorithm proceeds by:
        1. Representing each cluster as a PruningCluster dataclass
        2. Searching nearest-neighbor cluster pairs in embedding space
        3. Applying increasingly strict biological/statistical checks
        4. Merging the first failing pair per iteration
        5. Repeating until no merges remain

    This greedy strategy ensures that only locally similar clusters are tested.

    Parameters
    ----------
    membership : np.ndarray
        Initial cluster labels for each observation.
    gexpr : csr_matrix
        Sparse gene expression matrix (cells × genes), non log-transformed.
    embedding : np.ndarray
        Low-dimensional embedding used for geometric proximity search.
    strata : np.ndarray
        Stratification variable (e.g. sample, batch, donor) for Van Elteren test.
    nn_checks : int
        Number of nearest neighbor clusters to test per cluster.
    min_cluster_size : int
        Minimum cluster size before automatic merge.
    max_stratum_size : int
        Maximum number of cells sampled per stratum.
    min_stratum_size : int
        Minimum number of cells required per stratum.
    min_strata : int
        Minimum number of shared strata between clusters.
    min_pct_expression : float
        Minimum expression frequency for a gene to be considered.
    min_specificity : float
        Minimum combined fold-change and prevalence difference threshold.
    significance_threshold : float
        Family-wise p-value threshold (Bonferroni corrected internally).
    min_marker_genes : int
        Minimum number of significant genes required to keep clusters separate.
    symmetric: bool (default: True)
        Minimum marker genes is required of both positive and negative markers.

    Returns
    -------
    np.ndarray
        New cluster labels after pruning and merging.

    Notes
    -----
    • Uses Van Elteren stratified rank test for differential expression  
    • Applies Bonferroni correction across tested genes  
    • Merges are weighted by cluster size to preserve centroids and statistics  

    The algorithm terminates when all tested neighboring cluster pairs are
    statistically and biologically distinguishable.

    Complexity is roughly:

        O(K · nn_checks · DE_cost)

    where K is number of clusters.

    Practical performance is dominated by dense conversion for DE testing
    and multiseed leiden computation.
    """
    
    unique_clusters = np.unique(membership)
    
    # Setup ---------------------------------------------------------------------- #
    # Create PruningCluster dataclasses per cluster
    
    pruning_clusters = []
    for i, cluster in enumerate(unique_clusters):
        indices = np.where(membership == cluster)[0]
        centroid = embedding[indices].mean(axis=0)
        mean_gexpr = gexpr[indices].mean(axis=0).A1
        pct_gexpr = (gexpr[indices] > 0).mean(axis=0).A1
        
        pruning_clusters.append(
            PruningCluster(
                name=str(i),
                indices=indices,
                centroid=centroid,
                pct_gexpr=pct_gexpr,
                mean_gexpr=mean_gexpr,
            )
        )
    
    passed_pairs = set()
    any_merge_made = True
    
    while any_merge_made:
    
        centroids = np.array([c.centroid for c in pruning_clusters])
        pairs, _ = _get_prune_pairs(centroids, nn_checks=nn_checks)

        merge = False

        for pair in pairs:
            a, b = pair
            clusterA = pruning_clusters[a]
            clusterB = pruning_clusters[b]
            
            named_pair = frozenset([clusterA.name, clusterB.name])
            
            # No point recomputing already passed pairs
            if named_pair in passed_pairs:
                continue
            
            # Check 1: ------------------------------------------------------------------- #
            # Cluster Size
            
            if (clusterA.size < min_cluster_size) or (clusterB.size < min_cluster_size):
                merge = True
                print(f"{clusterA.name} or {clusterB.name} merged due to insufficient cluster size {clusterA.size}/{clusterB.size}.")
                break
            
            # Check 2: ------------------------------------------------------------------- #
            # Enough Marker Genes
            
            pct_A, pct_B = clusterA.pct_gexpr, clusterB.pct_gexpr
            pct_mask = (pct_A >= min_pct_expression) | (pct_B >= min_pct_expression)
            
            # Assumed expression is not log1p normalized
            log2_foldchange = np.log2(clusterA.mean_gexpr + (1/clusterA.size)) - np.log2(clusterB.mean_gexpr + (1/clusterB.size))
            
            rel_pct_diff = np.abs(pct_A - pct_B) / (pct_A + pct_B + 1e-10)
            specificity = rel_pct_diff * np.abs(log2_foldchange)
            
            specificity_mask = specificity >= min_specificity
            
            # Gene to be tested
            valid_genes = pct_mask & specificity_mask
            
            if sum(valid_genes) < min_marker_genes: # Preflight check
                merge = True
                print(f"{clusterA.name} and {clusterB.name} have insufficient valid genes.")
                break
            
            strata_A = strata[clusterA.indices]
            strata_B = strata[clusterB.indices]
            
            X_sub = gexpr[:, valid_genes]
            X_A = X_sub[clusterA.indices, :]
            X_B = X_sub[clusterB.indices, :]
            
            try:
                pvalues = van_elteren_test(
                    X=X_A, 
                    Y=X_B, 
                    X_strata=strata_A, 
                    Y_strata=strata_B,
                    min_stratum_size=min_stratum_size,
                    max_stratum_size=max_stratum_size,
                    min_strata=min_strata,
                    random_state=random_state
                )
            except ValueError: # insufficient overlapping strata between clusters
                merge = True
                print(f"{clusterA.name} and {clusterB.name} merged due to insufficient valid overlapping strata.")
                break
            
            # Bonferroni correction for multiple testing
            adj_significance_threshold = significance_threshold / np.sum(valid_genes)
            

            if symmetric:
                marker_sign_vec = np.sign(log2_foldchange[valid_genes])
                
                n_markers_A = np.sum(pvalues[marker_sign_vec == 1] < adj_significance_threshold)
                n_markers_B = np.sum(pvalues[marker_sign_vec == -1] < adj_significance_threshold)
                
                has_enough_A = n_markers_A >= min_marker_genes
                has_enough_B = n_markers_B >= min_marker_genes
                
                if not (has_enough_A and has_enough_B):
                    merge = True
                    print(f"{clusterA.name} vs {clusterB.name} FAILED with {n_markers_A} vs {n_markers_B} markers")
                    break
            else:
                
                n_markers = np.sum(pvalues < adj_significance_threshold)
                if n_markers < min_marker_genes:
                    merge = True
                    print(f"{clusterA.name} vs {clusterB.name} FAILED with {n_markers} markers")
                    break
            
            # If reached here then pairs are distinct (not pruned)
            passed_pairs.add(named_pair) # avoid recompute during loop
            if symmetric:
                print(f"{clusterA.name} vs {clusterB.name} passed with {n_markers_A} vs {n_markers_B} markers")
            else:
                print(f"{clusterA.name} vs {clusterB.name} passed with {n_markers} markers")            
            
        # On break, the failed pair is still `pair`
        if merge:
            i, j = sorted(pair, reverse=True) # pop order fix
            
            cluster_i = pruning_clusters.pop(i)
            cluster_j = pruning_clusters.pop(j)
            
            w_i, w_j = cluster_i.size, cluster_j.size
            w_t = w_i + w_j
            w_i, w_j = w_i / w_t, w_j / w_t
            
            cluster_ij = PruningCluster(
                name = "_".join([cluster_i.name, cluster_j.name]),
                indices = np.sort(np.union1d(cluster_i.indices, cluster_j.indices)),
                centroid = (w_i * cluster_i.centroid) + (w_j * cluster_j.centroid),
                pct_gexpr = (w_i * cluster_i.pct_gexpr) + (w_j * cluster_j.pct_gexpr),
                mean_gexpr = (w_i * cluster_i.mean_gexpr) + (w_j * cluster_j.mean_gexpr),
            )
            
            pruning_clusters.append(cluster_ij)
            
            print(f"\t>>{cluster_i.name} and {cluster_j.name} merged to make {cluster_ij.name}")
            
            merge = False
            
        else:
            any_merge_made = False
    
    # Rebuild membership vector from dataclasses --------------------------------- #
    
    pruned_consensus = np.full_like(membership, fill_value=-1)
    for i, cluster in enumerate(pruning_clusters):
        pruned_consensus[cluster.indices] = i

    return pruned_consensus
    

# ============================================================================ #
#                              Helper Function(s)                              #
# ============================================================================ #

# def _strata_subsample_indices(labels, indices, min_counts=10, max_counts=10_000, random_state=0):
#     """ [ChatGPT]
#     Subsample indices within each stratum while removing underrepresented strata.

#     Rare strata are discarded and abundant strata are randomly downsampled.
#     """
    
#     rng = np.random.default_rng(random_state)
    
#     subset_labels = labels[indices]
#     unique_labels, counts = np.unique(subset_labels, return_counts=True)
#     label_to_count = dict(zip(unique_labels, counts))
    
#     # Filter rare labels
#     mask = np.array([label_to_count[label] >= min_counts for label in subset_labels])
#     filtered_indices = indices[mask]
    
#     # Subsample frequent labels
#     final_indices = []
#     for label in np.unique(labels[filtered_indices]):
#         label_indices = filtered_indices[labels[filtered_indices] == label]
#         if len(label_indices) > max_counts:
#             final_indices.append(
#                 rng.choice(label_indices, size=max_counts, replace=False)
#             )
#         else:
#             final_indices.append(label_indices)
    
#     return np.concatenate(final_indices) if final_indices else np.array([], dtype=int)

def _get_prune_pairs(centroids, nn_checks=2):
    """ [ChatGPT]
    Return unique nearest-neighbor cluster index pairs ordered by distance.

    Uses kNN search in embedding space and removes symmetric duplicates.
    """
    
    n_centroids = centroids.shape[0]
    
    # n clusters may be smaller than nn_checks
    k_eff = min(nn_checks + 1, n_centroids)   # +1 because it includes self

    nbrs = NearestNeighbors(n_neighbors=k_eff, metric="euclidean")
    nbrs.fit(centroids)

    knn_dists, knn_indices = nbrs.kneighbors(centroids)

    # 1. Drop self if included as the first neighbor
    knn_dists, knn_indices = knn_dists[:, 1:], knn_indices[:, 1:]

    n_points, n_neighbors = knn_indices.shape

    # 2. Flatten arrays to make pairs
    self_idx = np.repeat(np.arange(n_points), n_neighbors)
    neighbor_idx = knn_indices.ravel()
    distances = knn_dists.ravel()

    # 3. Sort by distance first so smallest distance appears first
    sort_idx = np.argsort(distances)
    self_idx = self_idx[sort_idx]
    neighbor_idx = neighbor_idx[sort_idx]
    distances = distances[sort_idx]

    # 4. Make symmetric edges (min, max) to treat (a, b) == (b, a)
    pairs = np.stack([self_idx, neighbor_idx], axis=1)
    pairs = np.sort(pairs, axis=1)

    # 5. Deduplicate while preserving order using structured array trick
    dtype = np.dtype([('a', pairs.dtype), ('b', pairs.dtype)])
    structured_pairs = pairs.view(dtype)
    _, unique_idx = np.unique(structured_pairs, return_index=True)
    unique_idx.sort()  # preserve distance order

    # 6. Reformat pairs to tuples for later hashing
    pairs = pairs[unique_idx]
    pairs = [tuple(map(int, pair)) for pair in pairs]
    
    return pairs, distances[unique_idx]
