from typing import Iterable

import numpy as np

from numpy.typing import NDArray

from ..neighbors import perplexity_graph, knn_graph
from ..utils import csr_row_normalize

# def lisi(
#     probs_csr: csr_matrix,
#     labels: NDArray | pd.Series,
#     scaled=True
# ) -> NDArray:
    
#     # Compute label probabilities ------------------------------------------------ #

#     unique_obs_labels, obs_label_vec = np.unique(labels, return_inverse=True)
#     n_obs = probs_csr.shape[0]
#     n_batches = len(unique_obs_labels)

#     # Map columns (obs indices) to batch label indices
#     label_indices = obs_label_vec[probs_csr.indices]

#     # Define the shape for the sparse category sum matrix
#     matrix_shape = (n_obs, n_batches)

#     # Construct sparse matrix of label probability sums (neat trick: duplicate coordinates are summed)
#     label_probs_csr = csr_matrix((probs_csr.data.copy(), label_indices, probs_csr.indptr.copy()), shape=matrix_shape) # must copy or it breaks - space magic....
#     label_probs_csr.sum_duplicates() # IMPORTANT

#     # Compute LISI score --------------------------------------------------------- #
    
#     # Step 1: square probabilities
#     label_probs_csr.data **= 2

#     # Step 2: Compute Simpson Index (SI) for each observation base on neighbor probs giving Local Simpson Index (LSI)
#     lsi_vec = label_probs_csr.sum(axis=1).A1 # Local Simpson Index per obs

#     # Step 3: Transform to Inverse Simpson Index to get Local Inverse Simpson Index (LISI)
#     lisi_vec = 1 / lsi_vec  # Inverse Simpson Index transform

#     # Step 4 (optional): scale scores betwen 0 and 1
#     if scaled:
#         nnz_per_row = probs_csr.indptr[1:] - probs_csr.indptr[:-1]
#         max_nnz = np.max(nnz_per_row) # largest observed neighborhood, if knn then this will be k

#         # Scale by the maximum theoretical value min(number of categories, n_neighbors)
#         n_batches = len(unique_obs_labels)
#         max_lisi = min(max_nnz, n_batches)

#         lisi_vec = (lisi_vec - 1) / (max_lisi - 1)  # subtract one such that one batch has diversity of 0

#     return lisi_vec

def calculate_lisi_scores(
    embedding: NDArray,
    labels: NDArray | Iterable,
    perplexity: int = 30,
    scaled: bool =True,
    k_factor: int =3,
    **knn_kwargs
) -> NDArray:
    """Compute Local Inverse Simpson's Index (LISI) score for all points in a given embedding.

    Parameters
    ----------
    embedding : NDArray
        Embedding of each point as a 2D numpy array (n_points, n_dims).
    labels : NDArray | pd.Series
        A vector of batch levels associated with each point.
    perplexity : int, optional
        Effect number of nearest neighbors to consider per neighbourhood (perplexity), by default 30
    scaled : bool, optional
        Scale scores to the range 0-1 by dividing by theoretical maximum min(no. of neighbors total, no. of unique labels), by default True
    k_factor : int, optional
        Number of nearest neighbors to use for perplexity graph, by default 3, so 3 x 30 = 90 neighbors

    Returns
    -------
    NDArray
        A 1D array of LISI scores for each point.
    
    Notes
    -----
    The function uses a k-nearest neighbors approach to compute the conditional probabilities
    required for the LISI score. The LISI score is then calculated by computing the
    Local Inverse Simpson Index (LISI) for each observation based on the k-nearest neighbors and
    their corresponding batch labels.
    
    """    
    
    # Construct vector of integer label assignment for each obs
    unique_obs_labels, obs_label_inv = np.unique(labels, return_inverse=True)
    
    # Number of K nearest neighbors to find
    n_neighbors = k_factor * perplexity
    
    # Compute knn distance matrix (csr)
    D = knn_graph(
        embedding,
        k=n_neighbors + 1, # +1 as self is dropped later
        **knn_kwargs
    )
    D.eliminate_zeros() # remove self (assumes no duplicate points)
    
    # Construct perplexity connectivity matrix
    P = perplexity_graph(
        D, 
        perplexity=perplexity
    )
    
    # Normalize row-wise to derive probabilities from connectivities
    csr_row_normalize(P)
    
    # Construct binary adjacency for obs --> label mapping
    n_labels = len(unique_obs_labels)
    A = np.eye(n_labels, dtype=bool)[obs_label_inv]

    # Compute probabiliy each label for each obs (sum of all P(neighbors) with that label)
    P_labels = P @ A

    # Compute the Local Inverse Simpson Index (LISI) score for each observation
    P_labels **= 2  # Square the batch probabilities
    lsi_vec = P_labels.sum(axis=1)
    lisi_vec = 1 / lsi_vec  # Inverse Simpson index

    if scaled:
        # Scale by the maximum theoretical value min(number of categories, n neighbors)
        n_batches = len(unique_obs_labels)
        max_lisi = min(n_neighbors, n_batches)

        lisi_vec = (lisi_vec - 1) / (max_lisi - 1)  # subtract one such that one batch has diversity of 0

    return lisi_vec