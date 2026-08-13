import numpy as np
from numba import njit, types, prange
from scipy.sparse import csr_matrix

# ============================================================================ #
#                              Exposed Function(s)                             #
# ============================================================================ #

def umap_graph(
    D: csr_matrix,
    k: int | float,
    local_connectivity: float = 1.0,
    bandwidth: float = 1.0
) -> csr_matrix:
    """
    Compute a UMAP-style fuzzy simplicial set (sparse connectivity graph) 
    from a CSR distance matrix.

    Parameters
    ----------
    D : csr_matrix, shape (n_samples, n_samples)
        Sparse distance matrix.
    k : int | float
        Target number of neighbors for the fuzzy simplicial set.
    local_connectivity : float, default 1.0
        Minimum number of nearest neighbors to guarantee connectivity (ρ), can be fractional.
    bandwidth : float, default 1.0
        Scales the target sum of kernel weights.

    Returns
    -------
    csr_matrix
        Sparse matrix of UMAP fuzzy simplicial set weights.
    
    Notes
    -----
    - Distances equal to zero (typically self-distances) are ignored and kernel weight set to 0.
    - rho computation requires the distances to be sorted; fractional `local_connectivity` is linearly interpolated.
    - Binary search adjusts sigma to match the target sum of kernel values (similar to UMAP's smooth_knn_dist).
    - Kernel values for distances <= rho are set to 1.
    """
    connectivities = _umap_csr(
        distances=D.data,
        indptr=D.indptr,
        target_neighbors=k,
        local_connectivity=local_connectivity,
        bandwidth=bandwidth
    )

    C = csr_matrix((connectivities, D.indices.copy(), D.indptr.copy()), shape=D.shape)
    return C


def umap_knn(
    knn_dists: np.ndarray,
    k: int | float,
    local_connectivity: float = 1.0,
    bandwidth: float = 1.0
) -> np.ndarray:
    """
    Compute a UMAP-style fuzzy simplicial set from a dense k-nearest-neighbor distance matrix.

    Parameters
    ----------
    knn_dists : np.ndarray, shape (n_samples, n_neighbors)
        Distances to k-nearest neighbors.
    k : int | float
        Target number of neighbors for the fuzzy simplicial set.
    local_connectivity : float, default 1.0
        Minimum number of nearest neighbors to guarantee connectivity (ρ), can be fractional.
    bandwidth : float, default 1.0
        Scales the target sum of kernel weights.

    Returns
    -------
    np.ndarray, shape (n_samples, n_neighbors)
        UMAP fuzzy membership weights.
    
    Notes
    -----
    - Distances equal to zero (typically self-distances) are ignored and kernel weight set to 0.
    - rho computation requires the distances to be sorted; fractional `local_connectivity` is linearly interpolated.
    - Binary search adjusts sigma to match the target sum of kernel values (similar to UMAP's smooth_knn_dist).
    - Kernel values for distances <= rho are set to 1.
    """
    knn_conn = _umap_knn(
        knn_dists=knn_dists,
        target_neighbors=k,
        local_connectivity=local_connectivity,
        bandwidth=bandwidth        
    )
    
    return knn_conn

# ============================================================================ #
#                              Helper Function(s)                              #
# ============================================================================ #

@njit(parallel=True)
def _umap_csr(
    distances: np.ndarray,
    indptr: np.ndarray,
    target_neighbors: int | float,
    local_connectivity: float = 1.0,
    bandwidth: float = 1.0
) -> np.ndarray:
    """
    Compute UMAP fuzzy simplicial set for CSR distances.

    Parameters
    ----------
    distances : np.ndarray
        Nonzero distances in CSR format.
    indptr : np.ndarray
        CSR index pointer array.
    target_neighbors : int | float
        Target number of neighbors.
    local_connectivity : float
        Local connectivity parameter (ρ).
    bandwidth : float
        Bandwidth scaling factor.

    Returns
    -------
    np.ndarray
        Fuzzy connectivity weights for each nonzero entry.
    """
    connectivities = np.zeros_like(distances)
    n_samples = indptr.size - 1

    for i in prange(n_samples):
        start = indptr[i]
        end = indptr[i + 1]
        if end > start:
            connectivities[start:end] = _umap_kernel_transform(
                distances=distances[start:end], 
                target_neighbors=target_neighbors,
                local_connectivity=local_connectivity,
                bandwidth=bandwidth
            )

    return connectivities

@njit(parallel=True)
def _umap_knn(
    knn_dists: np.ndarray,
    target_neighbors: int | float,
    local_connectivity: float = 1.0,
    bandwidth: float = 1.0
) -> np.ndarray:
    """
    Compute UMAP fuzzy simplicial set for dense k-nearest-neighbor distances.

    Parameters
    ----------
    knn_dists : np.ndarray, shape (n_samples, n_neighbors)
        Distances to k-nearest neighbors.
    k : int | float
        Number of neighbors (ignored internally, included for consistency).
    target_neighbors : int | float
        Target number of neighbors for fuzzy membership.
    local_connectivity : float
        Local connectivity parameter (ρ).
    bandwidth : float
        Bandwidth scaling factor.

    Returns
    -------
    np.ndarray
        Fuzzy membership weights for each sample and neighbor.
    """
    n_samples, n_neighbors = knn_dists.shape
    knn_conn = np.zeros((n_samples, n_neighbors), dtype=np.float64)

    for i in prange(n_samples):
        knn_conn[i] = _umap_kernel_transform(
            distances=knn_dists[i], 
            target_neighbors=target_neighbors,
            local_connectivity=local_connectivity,
            bandwidth=bandwidth
        )

    return knn_conn


@njit(
    fastmath=True,
    locals={
        "num_neighbors": types.int64,
        "target_sum": types.float64,
        "sigma_min": types.float64,
        "sigma_max": types.float64,
        "sigma": types.float64,
        "kernel_sum": types.float64,
        "neighbor_idx": types.int64,
        "distance": types.float64,
        "kernel_value": types.float64,
        "i": types.int64,
        "rho_idx": types.int64,
        "rho": types.float64,
        "alpha": types.float64
    }
)
def _umap_kernel_transform(
    distances: np.ndarray,
    target_neighbors: int | float,
    local_connectivity: float = 1.0,
    bandwidth: float = 1.0,
    tol: float = 1e-5,
    max_iter: int = 64
) -> np.ndarray:
    """
    Compute UMAP-style fuzzy simplicial set weights for a single distance vector.

    Parameters
    ----------
    distances : np.ndarray, shape (n_neighbors,)
        Distances for a single sample. Zeros (self-distances) are masked.
    target_neighbors : int | float
        Target number of neighbors controlling the fuzzy kernel sum.
    local_connectivity : float, default 1.0
        Number of nearest neighbors guaranteed to be connected (ρ), can be fractional.
    bandwidth : float, default 1.0
        Scales the target sum of kernel weights.
    tol : float, default 1e-5
        Convergence tolerance for the binary search of sigma.
    max_iter : int, default 64
        Maximum iterations for the sigma search.

    Returns
    -------
    np.ndarray, shape (n_neighbors,)
        Fuzzy membership weights for each neighbor. Zeros correspond to masked distances.

    Notes
    -----
    - Distances equal to zero (typically self-distances) are ignored and kernel weight set to 0.
    - rho computation requires the distances to be sorted; fractional `local_connectivity` is linearly interpolated.
    - Binary search adjusts sigma to match the target sum of kernel values (similar to UMAP's smooth_knn_dist).
    - Kernel values for distances <= rho are set to 1.
    """
    num_neighbors = distances.shape[0]

    # --- Step 1: compute rho with fractional local_connectivity ---
    sorted_distances = np.sort(distances[distances > 0.0])
    if local_connectivity <= 0.0 or sorted_distances.size == 0:
        rho = 0.0
    else:
        f = local_connectivity
        i = int(np.floor(f))
        alpha = f - i
        if i >= sorted_distances.size - 1:
            rho = sorted_distances[-1]
        else:
            rho = (1.0 - alpha) * sorted_distances[i] + alpha * sorted_distances[i + 1]

    # --- Step 2: binary search for sigma ---
    target_sum = np.log2(target_neighbors) * bandwidth
    sigma_min = 1e-8
    sigma_max = 1e5
    sigma = 1.0

    kernel_values = np.empty(num_neighbors, dtype=np.float64)

    for _ in range(max_iter):
        kernel_sum = 0.0
        for neighbor_idx in range(num_neighbors):
            distance = distances[neighbor_idx]
            if distance == 0.0:
                kernel_value = 0.0
            elif distance <= rho:
                kernel_value = 1.0
            else:
                kernel_value = np.exp(-(distance - rho) / sigma)
            kernel_values[neighbor_idx] = kernel_value
            kernel_sum += kernel_value

        if np.abs(kernel_sum - target_sum) < tol:
            break

        if kernel_sum > target_sum:
            sigma_max = sigma
            sigma = 0.5 * (sigma + sigma_min)
        else:
            sigma_min = sigma
            sigma = 0.5 * (sigma + sigma_max)

    # --- Step 3: final kernel weights ---
    for neighbor_idx in range(num_neighbors):
        distance = distances[neighbor_idx]
        if distance == 0.0:
            kernel_values[neighbor_idx] = 0.0
        elif distance <= rho:
            kernel_values[neighbor_idx] = 1.0
        else:
            kernel_values[neighbor_idx] = np.exp(-(distance - rho) / sigma)

    return kernel_values
