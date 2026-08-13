import numpy as np
from numba import njit, prange, types
from scipy.sparse import csr_matrix

# ============================================================================ #
#                              Exposed Function(s)                             #
# ============================================================================ #

def perplexity_graph(
    D: csr_matrix, 
    perplexity: float
) -> csr_matrix:
    """
    Compute a Gaussian perplexity-weighted graph from a CSR distance matrix.

    Parameters
    ----------
    D : csr_matrix
        Sparse CSR distance matrix of shape (n_samples, n_samples).
    perplexity : float
        Target perplexity controlling effective neighborhood size.

    Returns
    -------
    csr_matrix
        Sparse CSR matrix of the same shape as D where each row contains
        Gaussian kernel weights normalized to match the target perplexity.
        Zero distances (self-distances) are ignored.
    """
    connectivities = _perplexity_csr(
        distances=D.data,
        indptr=D.indptr,
        target_perplexity=perplexity
    )

    C = csr_matrix((connectivities, D.indices.copy(), D.indptr.copy()), shape=D.shape)
    return C


def perplexity_knn(
    knn_dists: np.ndarray, 
    perplexity: float
) -> np.ndarray:
    """
    Compute a Gaussian perplexity-weighted connectivity matrix for dense KNN distances.

    Parameters
    ----------
    knn_dists : np.ndarray, shape (n_samples, n_neighbors)
        Pairwise distances for the k-nearest neighbors of each sample.
    perplexity : float
        Target perplexity controlling effective neighborhood size.

    Returns
    -------
    np.ndarray, shape (n_samples, n_neighbors)
        Gaussian kernel weights for each KNN row normalized to match the
        target perplexity.
    """
    knn_conn = _perplexity_knn(
        knn_dists=knn_dists,
        target_perplexity=perplexity
    )
    return knn_conn


# ============================================================================ #
#                              Helper Function(s)                              #
# ============================================================================ #

@njit(parallel=True)
def _perplexity_csr(
    distances: np.ndarray,
    indptr: np.ndarray,
    target_perplexity: float
) -> np.ndarray:
    """
    Compute row-wise Gaussian kernel weights for a CSR-style sparse distance matrix.

    Parameters
    ----------
    distances : np.ndarray
        Flattened non-zero distance values from a CSR matrix.
    indptr : np.ndarray
        CSR index pointer array.
    target_perplexity : float
        Target perplexity controlling effective neighborhood size.

    Returns
    -------
    np.ndarray
        Gaussian kernel weights for all non-zero distances.
    """
    connectivities = np.zeros_like(distances)
    n_samples = indptr.size - 1

    for i in prange(n_samples):
        start = indptr[i]
        end = indptr[i + 1]
        if end > start:  # skip empty rows
            connectivities[start:end] = _perplexity_transform(distances[start:end], target_perplexity)

    return connectivities


@njit(parallel=True)
def _perplexity_knn(
    knn_dists: np.ndarray, 
    target_perplexity: float
) -> np.ndarray:
    """
    Compute row-wise Gaussian kernel weights for dense KNN distances.

    Parameters
    ----------
    knn_dists : np.ndarray, shape (n_samples, n_neighbors)
        KNN distance matrix.
    target_perplexity : float
        Target perplexity controlling effective neighborhood size.

    Returns
    -------
    np.ndarray, shape (n_samples, n_neighbors)
        Gaussian kernel weights for each row.
    """
    n_samples, n_neighbors = knn_dists.shape
    knn_conn = np.zeros((n_samples, n_neighbors), dtype=np.float64)

    for i in prange(n_samples):
        knn_conn[i] = _perplexity_transform(knn_dists[i], target_perplexity)

    return knn_conn


@njit(
    fastmath=True,
    locals={
        "target_entropy": types.float64,
        "num_neighbors": types.int64,
        "sigma_min": types.float64,
        "sigma_max": types.float64,
        "sigma": types.float64,
        "inv_two_sigma2": types.float64,
        "kernel_sum": types.float64,
        "entropy": types.float64,
        "neighbor_idx": types.int64,
        "distance": types.float64,
        "kernel_value": types.float64,
        "prob": types.float64,
        "_": types.int64
    }
)
def _perplexity_transform(
    distances: np.ndarray,
    target_perplexity: float,
    tol: float = 1e-5,
    max_iter: int = 64
) -> np.ndarray:
    """
    Compute Gaussian kernel weights for a 1D distance vector that match a target perplexity.

    Parameters
    ----------
    distances : np.ndarray, shape (n_neighbors,)
        Distances for a single sample.
    target_perplexity : float
        Desired perplexity controlling effective number of neighbors.
    tol : float, optional
        Convergence tolerance for the binary search, by default 1e-5
    max_iter : int, optional
        Maximum number of iterations for the binary search, by default 60

    Returns
    -------
    np.ndarray, shape (n_neighbors,)
        Unnormalized Gaussian kernel weights corresponding to `distances`.
    """
    if distances.size == 0:
        return np.empty(0, dtype=np.float64)

    target_entropy = np.log(target_perplexity)
    num_neighbors = distances.shape[0]

    sigma_min = 1e-8
    sigma_max = 1e5
    sigma = 1.0

    kernel_values = np.empty(num_neighbors, dtype=np.float64)

    for _ in range(max_iter):
        inv_two_sigma2 = 1.0 / (2.0 * sigma * sigma)

        kernel_sum = 0.0
        for neighbor_idx in range(num_neighbors):
            distance = distances[neighbor_idx]
            if distance == 0.0:
                kernel_values[neighbor_idx] = 0.0
            else:
                kernel_value = np.exp(-distance * distance * inv_two_sigma2)
                kernel_values[neighbor_idx] = kernel_value
                kernel_sum += kernel_value

        entropy = 0.0
        if kernel_sum > 0.0:
            for neighbor_idx in range(num_neighbors):
                kernel_value = kernel_values[neighbor_idx]
                if kernel_value > 0.0:
                    prob = kernel_value / kernel_sum
                    entropy -= prob * np.log(prob)

        if np.abs(entropy - target_entropy) < tol:
            break

        if entropy > target_entropy:
            sigma_max = sigma
            sigma = 0.5 * (sigma + sigma_min)
        else:
            sigma_min = sigma
            sigma = 0.5 * (sigma + sigma_max)

    # final unnormalized Gaussian transform
    inv_two_sigma2 = 1.0 / (2.0 * sigma * sigma)
    for neighbor_idx in range(num_neighbors):
        distance = distances[neighbor_idx]
        if distance == 0.0:
            kernel_values[neighbor_idx] = 0.0
        else:
            kernel_values[neighbor_idx] = np.exp(-distance * distance * inv_two_sigma2)

    return kernel_values

