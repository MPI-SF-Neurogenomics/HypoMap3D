import numpy as np

from numba import njit
from scipy.sparse import csr_matrix

@njit
def separation_ratio_csr(
    data: np.ndarray,
    indices: np.ndarray,
    indptr: np.ndarray,
    labels: np.ndarray
) -> float:
    """
    Compute the separation ratio for a graph represented in CSR format.

    The separation ratio is defined as the ratio of the mean distance
    between points of different labels to the mean distance between
    points of the same label:

        separation_ratio = mean_between / mean_within

    Parameters
    ----------
    data : np.ndarray
        1D array of non-zero distances in the CSR matrix.
    indices : np.ndarray
        Column indices corresponding to `data` in CSR format.
    indptr : np.ndarray
        Row pointer indices for CSR format.
    labels : np.ndarray
        1D array of integer labels for each row/cell.

    Returns
    -------
    float
        Separation ratio between classes. Values >1 indicate that
        inter-class distances are larger than intra-class distances.
    """
    n = len(labels)
    
    sum_within = 0.0
    count_within = 0
    sum_between = 0.0
    count_between = 0
    
    for i in range(n):
        row_start = indptr[i]
        row_end = indptr[i+1]
        
        for idx in range(row_start, row_end):
            j = indices[idx]
            if i == j:
                continue  # skip self-distance
            
            d = data[idx]
            if labels[i] == labels[j]:
                sum_within += d
                count_within += 1
            else:
                sum_between += d
                count_between += 1
    
    mean_within = sum_within / count_within
    mean_between = sum_between / count_between
    
    return mean_between / mean_within


def compute_neighborhood_jaccard(
    A: csr_matrix,
    B: csr_matrix
) -> np.ndarray:
    """
    Compute row-wise Jaccard similarity between two CSR sparse matrices.

    For each row (cell), the Jaccard similarity is:

        J = |neighbors_A ∩ neighbors_B| / |neighbors_A ∪ neighbors_B|

    Parameters
    ----------
    A : csr_matrix
        First sparse adjacency or KNN matrix.
    B : csr_matrix
        Second sparse adjacency or KNN matrix.

    Returns
    -------
    np.ndarray
        1D array of Jaccard similarity values, one per row. Values
        range from 0 (no shared neighbors) to 1 (identical neighbors).
    """
    # Convert to binary adjacency
    A_bin = (A > 0).astype(int)
    B_bin = (B > 0).astype(int)

    # Compute intersection per row
    intersection = A_bin.multiply(B_bin).sum(axis=1).A1  # .A1 converts (n,1) to 1D array

    # Compute union per row
    union = A_bin.sum(axis=1).A1 + B_bin.sum(axis=1).A1 - intersection

    # Compute Jaccard similarity, safely avoiding division by zero
    jaccard = np.zeros_like(intersection, dtype=float)
    nonzero_mask = union != 0
    jaccard[nonzero_mask] = intersection[nonzero_mask] / union[nonzero_mask]

    return jaccard

def neighborhood_agreement(
    D: csr_matrix,
    memberships: np.ndarray
):
    """ [ChatGPT]
    Compute the fraction of nearest neighbors sharing the same cluster for each node.

    Given a k-nearest neighbor (kNN) graph in sparse CSR format and an array of cluster 
    memberships, this function calculates, for each node, the proportion of its neighbors 
    that belong to the same cluster. This is commonly used as a local cluster agreement 
    metric in clustering evaluation.

    Parameters
    ----------
    D : scipy.sparse.csr_matrix, shape (n_nodes, n_nodes)
        Sparse adjacency matrix of the kNN graph. Nonzero entries indicate edges (neighbors) 
        between nodes. Typically, the diagonal should be zero (no self-loops).
    memberships : np.ndarray, shape (n_nodes,)
        Array of cluster labels for each node. Labels are automatically remapped to 
        contiguous integers internally.

    Returns
    -------
    pct_agreement : np.ndarray, shape (n_nodes,)
        Fraction of neighbors that share the same cluster as the node itself. Each entry
        is in [0, 1], with 1 indicating perfect agreement (all neighbors in the same cluster).

    Notes
    -----
    - The function assumes each row of `D` represents the neighbors of a node and that 
      the graph is a kNN graph with roughly equal numbers of neighbors per node.
    - `D` must be in CSR format for efficient row-wise access.
    - Self-neighbors (defined by distance of zero) are ignored.
    """
    n_obs, _ = D.shape
    
    _, memberships = np.unique(memberships, return_inverse=True)
    
    self_mask = D.data == 0
    knn_indices = D.indices[~self_mask].reshape((n_obs, -1))
    
    knn_memberships = memberships[knn_indices]
    
    match_mask = knn_memberships == memberships[:, None]
    pct_agreement = match_mask.mean(axis=1)
    
    return pct_agreement