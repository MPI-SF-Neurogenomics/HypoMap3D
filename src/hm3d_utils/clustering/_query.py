import logging 
import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
import torch
from typing import List, Union

from sklearn.preprocessing import OrdinalEncoder, OneHotEncoder

from ..neighbors import knn_graph, umap_graph
from ..utils import csr_row_normalize

use_gpu = torch.cuda.is_available()

if use_gpu:
    try:
        from cuml.neighbors import NearestNeighbors as CuNearestNeighbors
    except ImportError:
        print("CUDA is available but cuML is not installed. Falling back to CPU.")
        from sklearn.neighbors import NearestNeighbors as SkNearestNeighbors
        use_gpu = False
else:
    from sklearn.neighbors import NearestNeighbors as SkNearestNeighbors

def transfer_labels(
    X_ref: np.ndarray,
    X_query: np.ndarray,
    lbls_ref: np.ndarray,
    n_neighbors: int = 50,
    use_gpu: bool = True,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Transfer categorical labels from a reference dataset to a query dataset using 
    nearest-neighbor propagation.

    This function computes a k-nearest neighbors (KNN) graph between reference 
    points and query points, transforms neighbor distances into connectivity weights, 
    and propagates one-hot encoded labels from the reference to the query points. 
    The most probable label is assigned to each query point.

    Parameters
    ----------
    X_ref : np.ndarray of shape (n_ref_samples, n_features)
        Feature matrix of the reference dataset.
    X_query : np.ndarray of shape (n_query_samples, n_features)
        Feature matrix of the query dataset to which labels will be transferred.
    lbls_ref : np.ndarray of shape (n_ref_samples,)
        Array of categorical labels corresponding to `X_ref`.
    n_neighbors : int, optional, default=50
        Number of nearest neighbors to consider in the KNN graph.
    use_gpu : bool, optional, default=True
        Whether to use GPU acceleration for nearest-neighbor computation if available.

    Returns
    -------
    lbls_query : np.ndarray of shape (n_query_samples,)
        Predicted labels for the query dataset, inferred from the reference labels.
    lbls_probs : np.ndarray of shape (n_query_samples,)
        Probability of the assigned label for each query point, based on 
        normalized neighbor connectivity.

    Notes
    -----
    - Reference labels are encoded with `OrdinalEncoder` for integer representation 
      and then transformed into one-hot matrices using `OneHotEncoder`.
    - Distances between query and reference points are converted into weights 
      via `hm3d_utils.nn.umap_graph` and normalized row-wise to form probabilities.
    - Label probabilities for query points are computed as the weighted sum of 
      reference labels.
    - The most probable label is selected as the final prediction.
    
    Examples
    --------
    >>> lbls_query, lbls_probs = transfer_labels(
    ...     X_ref=ref_features, 
    ...     X_query=query_features, 
    ...     lbls_ref=ref_labels,
    ...     n_neighbors=30
    ... )
    """
    # Convert to numpy array for onehot encoder (e.g. handle DataFrame)
    lbls_ref = np.asarray(lbls_ref)[:, None]

    # Encoder label categories as integer labels
    lbl_enc = OrdinalEncoder(dtype=int)
    lbl_ref_encoded = lbl_enc.fit_transform(lbls_ref)

    # Construct KNN graph of nearest reference neighbor to each query point
    D = knn_graph( # (n_query x n_ref)
        X_ref,
        k=n_neighbors,
        query=X_query,
        use_gpu=use_gpu
    )

    # Transform distances to connectivity weights
    C_nn = umap_graph( # (n_query x n_ref)
        D,
        k=n_neighbors, # target k
    )

    # Convert weights to probabiltiies
    P_nn = C_nn.copy()
    csr_row_normalize(P_nn) # (n_query x n_ref)

    # Setup encoder for conversion of integer category labels to sparse adjacency matrices
    oh_enc = OneHotEncoder(dtype=bool)
    oh_enc.fit(lbl_ref_encoded)
    A_lbls_ref = oh_enc.transform(lbl_ref_encoded)

    P_lbls_ref = P_nn @ A_lbls_ref

    max_lbls = P_lbls_ref.argmax(axis=1).A1
    max_path_probs = P_lbls_ref.max(axis=1).todense().A1

    lbls_query = lbl_enc.inverse_transform(max_lbls[:, None]).flatten()
    lbls_probs = max_path_probs
    
    return lbls_query, lbls_probs

def PropagateLabels(reference: ad.AnnData, query: ad.AnnData, label: str, 
        embed: str ='X_embed', n_neighbors:int=50,):
    """ Wrapper for new label transfer function transfer_labels"""

    query_label, label_confidence = transfer_labels(
            X_ref=reference.obsm[embed], 
            X_query=query.obsm[embed], 
            lbls_ref=reference.obs[label].values, 
            n_neighbors=n_neighbors, use_gpu=use_gpu)

    query_label = pd.Series(query_label, index=query.obs.index, dtype="category")
    label_confidence = pd.Series(label_confidence, index=query.obs.index, dtype='float32')
    
    return query_label, label_confidence


def propagate_labels_within_dataset(adata: ad.AnnData, propagation_label: str, 
                                    propagation_key: str | None = None,) -> None:
    """
    Propagates missing labels within an AnnData object using a label propagation method.

    This function identifies unlabelled entries in the `adata.obs` DataFrame under the specified 
    `propagation_label` column and attempts to infer those labels based on the labelled entries. 
    The inferred labels and their associated confidence scores are stored back in `adata.obs`.

    Parameters:
    ----------
    adata : ad.AnnData
        The annotated data matrix containing both labelled and unlabelled data.
    propagation_label : str
        The column name in `adata.obs` indicating which label to propagate.
    propagation_key : str | None, optional
        A base key used to name the output columns for propagated labels and their confidence.
        If None, defaults to the name of the `propagation_label`.

    Returns:
    -------
    None
        The function modifies `adata.obs` in-place by adding two new columns:
        - '{propagation_key}_propagated': containing the original and inferred labels.
        - '{propagation_key}_confidence': containing the confidence scores for the propagated labels.
    """

    if propagation_key == None:
        propagation_key = propagation_label

    mask = adata.obs[propagation_label].isna()

    unlabelled_view = adata[mask]
    labelled_view = adata[~mask]

    propagated_labels, propagated_confidence = transfer_labels(X_ref=labelled_view.obsm["X_embed"], 
                                X_query=unlabelled_view.obsm["X_embed"], 
                                lbls_ref=labelled_view.obs[propagation_label].values, 
                                n_neighbors=50, use_gpu=use_gpu)

    adata.obs[f'{propagation_key}_propagated'] = adata.obs[propagation_label]
    adata.obs[f'{propagation_key}_confidence'] = pd.Series(None, dtype=np.float32)

    adata.obs.loc[unlabelled_view.obs.index, f'{propagation_key}_propagated'] = propagated_labels
    adata.obs.loc[unlabelled_view.obs.index, f'{propagation_key}_confidence'] = propagated_confidence


def project_query(
    ref_adata: ad.AnnData,
    qry_adata: ad.AnnData,
    path_keys: list | tuple,
    k: int = 50,
    use_rep: str = 'X_embed',
    use_gpu: bool = True
) -> tuple[np.ndarray, sp.csr_matrix]:
    """
    Project hierarchical path annotations from a reference AnnData object onto
    query cells using a weighted bipartite shared-nearest-neighbor (SNN) graph
    constructed in a shared embedding space.

    The method operates in four stages:

    1. Construction of reference and query k-nearest-neighbor graphs
    2. Computation of bipartite shared-neighbor similarities
    3. Conversion of similarities into transition probabilities
    4. Probabilistic projection of hierarchical paths onto query cells

    Notes
    -----
    - The method assumes that both reference and query embeddings reside in
      the same latent space.
    - Path probabilities are heuristic consensus scores derived from
      neighborhood composition and should not be interpreted as conditional
      probabilities.

    Parameters
    ----------
    ref_adata : ad.AnnData
        Reference AnnData object containing:

        - ``.obsm[use_rep]``: embedding representation
        - ``.obs[path_keys]``: ordered hierarchical labels
        - optionally ``.obsp[f'{k}nn_adjacencies']`` containing a precomputed
          reference kNN adjacency matrix

    qry_adata : ad.AnnData
        Query AnnData object containing ``.obsm[use_rep]``.

    path_keys : list | tuple
        Ordered column names defining the hierarchical annotation structure.

        Example
        -------
        ``['Lineage', 'Class', 'Subclass', 'State']``

    k : int, optional
        Number of nearest neighbors used during graph construction,
        by default ``50``.

    use_rep : str, optional
        Key in ``.obsm`` specifying the shared embedding representation,
        by default ``'X_embed'``.

    use_gpu : bool, optional
        Whether to use GPU-accelerated nearest-neighbor search via
        ``cuml.neighbors`` when available, by default ``True``.

    Returns
    -------
    pathmat : np.ndarray
        Array of shape ``(n_paths, n_levels)`` containing all unique
        hierarchical paths observed in the reference dataset.

    qry_P_path : scipy.sparse.csr_matrix
        Sparse matrix of shape ``(n_query, n_paths)`` containing normalized
        path probabilities for each query cell.
    """
    
    ref_X = ref_adata.obsm[use_rep]
    ref_labelmat = ref_adata.obs[path_keys]
    
    qry_X = qry_adata.obsm[use_rep]
    
    # ---------------------------- Reference Neighbors --------------------------- #
    
    # In order to compute bipartite SNN we need reference neighbourhoods
    if f'{k}nn_adjacencies' not in ref_adata.obsp.keys():
        logging.warning(
            f"""
            Precomputed reference neighbor adjacencies, `.obsp['{k}nn_adjacencies']`, not found. 
            Computing reference {k} neighbors and adding to `.obsp['{k}nn_adjacencies']`.
            """
        )
        
        D = knn_graph(
            ref_adata.obsm['X_embed'],
            k=k+1, # +1 for self reference
            use_gpu=use_gpu
        )
        
        ref_adata.obsp[f'{k}nn_adjacencies'] = (D > 0) # convert to boolean and drop self (d=0)
    
    ref_A = ref_adata.obsp[f'{k}nn_adjacencies']
    
    # ------------------------------ Query Neighbors ----------------------------- #
    
    # Compute indices of nearest reference points for all query points as an adjacency matrix (qry_A)
    D = knn_graph(
        ref_X,
        k=k,
        query=qry_X,
        use_gpu=use_gpu
    )
    qry_A = (D > 0).astype(np.int32) # unfortunately can't both be boolean, qry_A is smaller so it is made integer
    
    # ----------------------- Weighted KNN (bipartite SNN) ----------------------- #
    
    # Compute number of shared nearest neighbors between qry and ref neighborhoods
    qry_P_nn = qry_A @ ref_A.T

    # Mask out non-neighbors with connection after sparse multiplication
    qry_P_nn = qry_P_nn.multiply(qry_A)

    # Compute jaccard - union = 100 + 100 - intersection
    qry_P_nn.data = (qry_P_nn.data / (2*k - qry_P_nn.data)).astype(np.float32)

    # Normalize to probabilities
    csr_row_normalize(qry_P_nn)
    
    # ----------------------------- Path Probabilties ---------------------------- #
    
    path_length = len(path_keys)

    # Encoder label categories as integer labels
    lbl_enc = OrdinalEncoder(dtype=int)
    ref_labelmat_encoded = lbl_enc.fit_transform(ref_labelmat)

    # Dynamically identify all possible paths in reference hierarchy
    paths_encoded = np.unique(ref_labelmat_encoded, axis=0)

    # Setup encoder for conversion of integer category labels to sparse adjacency matrices
    oh_enc = OneHotEncoder(dtype=bool)
    oh_enc.fit(paths_encoded)

    # Create adjacency matrices for membership of labels to each reference points and to each possible path
    ref_A_lbls = oh_enc.transform(ref_labelmat_encoded) # (n_ref x n_clusters)
    A_paths = oh_enc.transform(paths_encoded) # (n_paths x n_clusters)
    
    # Compute probability of all labels for each query point (majority vote structure)
    qry_P_lbls = qry_P_nn @ ref_A_lbls # (n_query x n_clusters)

    # Setup a binary variant for nonzero values
    qry_P_lbls_bin = (qry_P_lbls != 0).copy() # (n_query x n_clusters)

    # Convert to logspace for multiplicative probability along each path via summation
    qry_P_lbls.data = np.log(
        np.clip(qry_P_lbls.data, 1e-9, 0.99999) # P(x) = 1 breaks sparsity
    )
    qry_P_path = qry_P_lbls @ A_paths.T # (n_query x n_paths)
    
    # Exponentiate to escape to log space
    qry_P_path.data = np.exp(qry_P_path.data / path_length)

    # Any path that contianed a zero must be set to 0 probabiltiy (sparse malarky)
    coverage = qry_P_lbls_bin.astype(np.int32) @ A_paths.T
    valid = coverage == path_length # (n_query x n_paths)
    qry_P_path = qry_P_path.multiply(valid)

    # Normalize probabilty/energy across all paths
    csr_row_normalize(qry_P_path)

    # Convert encoded paths to original labels for return reference
    pathmat = lbl_enc.inverse_transform(paths_encoded)

    return pathmat, qry_P_path