from scipy.sparse import csr_matrix

from typing import Literal
from numpy.typing import NDArray

def knn_graph(
    reference: NDArray,
    k: int,
    query: NDArray = None,
    use_gpu: bool = True,
    **kwargs
) -> csr_matrix:
    """ [ChatGPT]
    Construct a k-nearest neighbors (kNN) graph as a sparse adjacency matrix.

    This function computes the kNN graph of a set of reference points, optionally
    including a separate set of query points. The graph can be computed using either
    scikit-learn or cuML depending on whether CUDA is enabled in `scutils.config`.

    Parameters
    ----------
    reference : NDArray
        A 1D or 2D array of reference points (embeddings/features).
    k : int
        Number of nearest neighbors to consider for each point.
    query : NDArray, optional
        Optional query points. If provided, the graph is computed between
        query points and reference points (query × reference). If None,
        the graph is computed among reference points (reference × reference).
    **kwargs
        Additional keyword arguments passed to the underlying kNN estimator
        (`cuml.neighbors.NearestNeighbors` if CUDA is enabled, otherwise
        `sklearn.neighbors.NearestNeighbors`).

    Returns
    -------
    csr_matrix
        Sparse adjacency matrix of shape (n_query, n_reference) or
        (n_reference, n_reference) if `query` is None. The matrix entries
        contain the distances to the nearest neighbors.
    
    Notes
    -----
    - Distances are used as edge weights in the adjacency matrix.
    - The underlying kNN algorithm depends on the availability of GPU support.
    """
    
    model = _knn_estimator(k, use_gpu, **kwargs)
    model.fit(reference)

    X = reference if query is None else query
    
    D = model.kneighbors_graph(X, mode='distance')
    D = _to_cpu(D)
    
    return D

def knn(
    reference: NDArray,
    k: int,
    query: NDArray = None,
    use_gpu: bool = True,
    **kwargs
) -> tuple[NDArray, NDArray]:
    """ [ChatGPT]
    Compute the k-nearest neighbors (kNN) of a set of points and return distances and indices.

    This function computes the k-nearest neighbors of reference points, optionally
    including a separate set of query points. The computation defaults to scikit-learn 
    or cuML depending on whether CUDA is enabled in `scutils.config`.

    Parameters
    ----------
    reference : NDArray
        A 1D or 2D array of reference points (embeddings/features).
    k : int
        Number of nearest neighbors to compute for each point.
    query : NDArray, optional
        Optional query points. If provided, neighbors are computed between
        query points and reference points (query × reference). If None,
        neighbors are computed among reference points (reference × reference).
    **kwargs
        Additional keyword arguments passed to the underlying kNN estimator
        (`cuml.neighbors.NearestNeighbors` if CUDA is enabled, otherwise
        `sklearn.neighbors.NearestNeighbors`).

    Returns
    -------
    tuple[NDArray, NDArray]
        - `distances`: Array of shape (n_query, k) containing the distances to the
          k nearest neighbors.
        - `indices`: Array of shape (n_query, k) containing the indices of the
          nearest neighbors in the reference set.
    
    Notes
    -----
    - Distances and indices are returned separately, unlike `knn_graph`, which returns
      a sparse adjacency matrix.
    - The underlying kNN algorithm depends on GPU availability.
    """
    
    model = _knn_estimator(k, use_gpu, **kwargs)
    model.fit(reference)
    
    X = reference if query is None else query
    
    knn_dist, knn_indices = model.kneighbors(X, return_distance=True)
    knn_dist = _to_cpu(knn_dist)
    knn_indices = _to_cpu(knn_indices)

    return knn_dist, knn_indices
    
# def perplexity_graph(
#     reference: NDArray,
#     p: int ,
#     k: int = None,
#     query: NDArray | None = None,
#     return_mtx: bool = False,
#     **knn_kwargs
# ) -> csr_matrix | tuple[NDArray, NDArray]:
#     """
#     Notes:
#     Returns a directed CSR matrix representaiton of perplexity graph
#     To convert to undirected use `scutils.nn.symmetrize_graph(P, 'union')` - MUST be ref x ref shape
#     """
#     k = k if k else 3 * p
#     knn_dists, knn_idx = knn_graph(
#         reference=reference, 
#         k=k, 
#         query=query, 
#         return_mtx=False,
#         **knn_kwargs
#         )
    
#     knn_dists = knn_dists.astype(np.float32, copy=False) # type requirement of cython perplexity search
#     P_cond = _knn_bs_perplexity(knn_dists, p)
    
#     if return_mtx:
#         shape = (reference.shape[0], reference.shape[0]) if query is None else (query.shape[0], reference.shape[0])
#         P = _knn_values_to_csr_matrix(P_cond, knn_idx, shape=shape)
#         P.eliminate_zeros()
#         return P
#     else:
#         return P_cond, knn_idx

# def snn_graph(
#     reference: NDArray,
#     k: int,
#     query: NDArray | None = None,
#     return_mtx: bool = False,
#     **knn_kwargs
# ) -> csr_matrix | tuple[NDArray, NDArray]:
#     """
#     Construct a shared nearest neighbor (SNN) graph.

#     The graph is built from a k-nearest neighbors (kNN) search, where edges are
#     weighted based on the number of shared neighbors (Jaccard similarity).  

#     Parameters
#     ----------
#     reference : NDArray
#         A 1D or 2D embedding of reference points.
#     k : int
#         Number of nearest neighbors to consider.
#     query : NDArray, optional
#         Optional query points. If provided, the graph is computed between
#         query x reference. If None (default), the graph is computed as
#         reference × reference.
#     return_mtx : bool, default=False
#         Determines the return format:
#         - If False, returns a tuple ``(knn_dists, knn_indices)``.
#         - If True, returns the adjacency matrix as a ``csr_matrix``.
#     **knn_kwargs
#         Additional keyword arguments passed to the underlying kNN model
#         (``cuml.neighbors.NearestNeighbors`` if CUDA is enabled,
#         otherwise ``sklearn.neighbors.NearestNeighbors``).

#     Returns
#     -------
#     csr_matrix or tuple of (NDArray, NDArray)
#         - If ``return_mtx=True``: a sparse adjacency matrix of shape
#           (ref × ref) or (query × ref).
#         - If ``return_mtx=False``: a tuple ``(knn_dists, knn_indices)``
#           containing neighbor distances and indices.

#     Notes
#     -----
#     - Recommended to use a larger ``k`` than typical UMAP defaults
#       (e.g., ``k=20`` instead of ``k=15``).
#     - The output graph is directed.
#     """

#     _, knn_idxs = knn_graph(
#         reference=reference, 
#         k=k,
#         query=query,
#         return_mtx=False,
#         **knn_kwargs
#         )
    
#     snn_dists = _calculate_shared_neighbors(knn_idxs)
    
#     if return_mtx:
#         shape = (query.shape[0], reference.shape[0]) if query else (reference.shape[0], reference.shape[0])
#         return _knn_values_to_csr_matrix(snn_dists, knn_idxs, shape=shape)
#     else:
#         return snn_dists, knn_idxs

def symmetrize_graph(
    graph_mtx: csr_matrix,
    type: Literal['mutual', 'union', 'avg', 'max', 'min'] = 'max',
) -> csr_matrix:
    """
    Symmetrize a directed graph adjacency matrix.

    Applies a chosen symmetrization strategy to a directed graph to produce
    an undirected graph representation.

    Parameters
    ----------
    graph_mtx : csr_matrix
        The input directed graph in CSR (Compressed Sparse Row) matrix format.
        Must be square.
    type : {'mutual', 'union', 'avg', 'max', 'min'}, default='max'
        The symmetrization method:
        - 'mutual' : Keep only edges that appear in both directions.
        - 'union'  : Keep edges that appear in either direction.
        - 'avg'    : Use the average of forward and reverse edge weights.
        - 'max'    : Use the maximum of forward and reverse edge weights.
        - 'min'    : Use the minimum of forward and reverse edge weights.

    Returns
    -------
    csr_matrix
        A symmetric adjacency matrix in CSR format representing the undirected graph.

    Raises
    ------
    ValueError
        If ``graph_mtx`` is not square.
    ValueError
        If an invalid ``type`` is provided.
    """

    if graph_mtx.shape[0] != graph_mtx.shape[1]:
        raise ValueError("Matrix must be square to symmetrize.")

    A = graph_mtx
    AT = A.T

    if type == 'mutual':
        sym = A.minimum(AT)
    elif type == 'union':
        sym = A + AT
    elif type == 'avg':
        sym = (A + AT) * 0.5
    elif type == 'max':
        sym = A.maximum(AT)
    elif type == 'min':
        sym = A.minimum(AT)
    else:
        raise ValueError(f"Unknown symmetrization type: {type}")

    return sym.to_csr()

# ## Helper Functions
# # --------------------------------------------------------------------------

# def _knn_values_to_csr_matrix(knn_dists: NDArray, knn_indices: NDArray, shape: tuple[int, int]):
#     """Convert knn_dists + knn_indices arrays into a CSR matrix format"""
#     n_samples, n_neighbors = knn_dists.shape

#     # Flatten the data and indices
#     data = knn_dists.ravel()
#     indices = knn_indices.ravel()

#     # Construct indptr
#     indptr = np.arange(0, n_samples * n_neighbors + 1, n_neighbors)

#     # Build the CSR matrix
#     graph_mtx = csr_matrix((data, indices, indptr), shape=shape)
#     graph_mtx.eliminate_zeros() # remove 0 data points - not automatic in construction
    
#     return graph_mtx

def _knn_estimator(k: int, use_gpu=False, **kwargs):
    """A helper function to use either a NearestNeigbbors model from sklearn or cuML dependant on GPU availability."""
    if use_gpu: # GPU neighborhood with cuML
        from cuml.neighbors import NearestNeighbors as cuNearestNeighbors
        model = cuNearestNeighbors(n_neighbors=k, **kwargs)
    else: # CPU neighborhood with sklearn
        from sklearn.neighbors import NearestNeighbors as skNearestNeighbors
        model = skNearestNeighbors(n_neighbors=k, **kwargs)

    return model

def _to_cpu(x):
    """Convert CuPy objects to CPU if needed."""
    return x.get() if hasattr(x, "get") else x

# @njit(parallel=True)
# def _calculate_shared_neighbors(
#     indices: NDArray,
#     prune_snn: float | None = 1/15,
# ) -> NDArray:
#     """
#     Compute shared nearest neighbor (SNN) weights using Jaccard similarity.

#     For each sample, this function measures the overlap between its k-nearest
#     neighbors and those of each neighbor. Edge weights are given by the Jaccard
#     index of the neighbor sets, with optional pruning.

#     Parameters
#     ----------
#     indices : NDArray of shape (n_samples, k)
#         Indices of the k-nearest neighbors for each sample.
#     prune_snn : float or None, default=1/15
#         Pruning threshold for edge weights:
#         - If not None, edges with Jaccard weight below this threshold are
#           set to 0 (pruned).
#         - If None, no pruning is applied.
#         By default, ``1/15`` is used (following Seurat’s recommendation).

#     Returns
#     -------
#     snn_values : NDArray of shape (n_samples, k)
#         Array of SNN weights corresponding to neighbor edges. Entries are
#         either the Jaccard weight or 0.0 if pruned.

#     Notes
#     -----
#     - The Jaccard weight is defined as:
#       ``overlap_size / (2 * k - overlap_size)``  
#       where ``overlap_size`` is the number of shared neighbors.
#     - Parallelized with Numba for efficiency.
#     """
#     n_samples, k = indices.shape
#     snn_values = np.zeros((n_samples, k), dtype=np.float64)

#     for i in prange(n_samples):
#         Ni = indices[i]
#         set_Ni = set(Ni)  # for O(1) lookups
#         for idx_j in range(k):
#             j = Ni[idx_j]
#             # count shared neighbours
#             overlap = 0
#             for x in indices[j]:
#                 if x in set_Ni:
#                     overlap += 1
#             # Jaccard
#             union_size = 2*k - overlap
#             w = overlap / union_size if union_size > 0 else 0.0
#             # Prune if enabled (defualt 1/15 same as Seruat)
#             if prune_snn is not None and w < prune_snn:
#                 snn_values[i, idx_j] = 0.0
#             else:
#                 snn_values[i, idx_j] = w

#     return snn_values
    
# Source - https://stackoverflow.com/a/53013222
# Posted by user5751
# Retrieved 2026-02-09, License - CC BY-SA 4.0

# def pairwise_jaccard_sparse_weighted(csr, epsilon, weight):
#     csr = scipy.sparse.csr_matrix(csr).astype(bool).astype(int)
#     csr_w = csr * scipy.sparse.diags(weight)

#     csr_rowsum = numpy.array(csr_w.sum(axis = 1)).flatten()
#     intrsct = csr.dot(csr_w.T)

#     rowsum_i = numpy.repeat(csr_rowsum, intrsct.getnnz(axis = 1))
#     unions = rowsum_i + csr_rowsum[intrsct.indices] - intrsct.data
#     dists = 1.0 - 1.0 * intrsct.data / unions

#     mask = (dists > 0) & (dists <= epsilon)
#     data = dists[mask]
#     indices = intrsct.indices[mask]

#     rownnz = numpy.add.reduceat(mask, intrsct.indptr[:-1])
#     indptr = numpy.r_[0, numpy.cumsum(rownnz)]

#     out = scipy.sparse.csr_matrix((data, indices, indptr), intrsct.shape)
#     return out
