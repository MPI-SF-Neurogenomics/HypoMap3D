import numpy as np

def grouped_means(data: np.ndarray, groups: np.ndarray, return_mapped: bool = False) -> np.ndarray:
    """
    Compute the mean of each cluster for 1D or ND data.

    Parameters
    ----------
    data : np.ndarray
        Array of shape (N,) for 1D data or (N, D) for multi-dimensional data.
    groups : np.ndarray
        Array of shape (N,), cluster IDs (integers or arbitrary labels).
    return_mapped : bool, default=False
        If True, return the per-point cluster mean mapped back to original points.
        If False, return cluster means only.

    Returns
    -------
    np.ndarray
        If return_mapped=False:
            Cluster means:
                - shape (K,) for 1D data
                - shape (K, D) for ND data
        If return_mapped=True:
            Array of same shape as `data`, with each point replaced by its cluster mean.
    """
    data = np.asarray(data)
    group_labels = np.asarray(groups)

    # Handle 1D data uniformly
    is_1d = data.ndim == 1
    if is_1d:
        data = data[:, None]  # make it 2D for consistent handling

    # Normalize labels to integers 0..K-1
    unique_labels, inverse_labels = np.unique(group_labels, return_inverse=True)
    K = len(unique_labels)

    # Allocate sums and counts
    cluster_sums = np.zeros((K, data.shape[1]), dtype=data.dtype)
    cluster_counts = np.zeros(K, dtype=np.int64)

    # Sum and count points per cluster
    np.add.at(cluster_sums, inverse_labels, data)
    np.add.at(cluster_counts, inverse_labels, 1)

    # Compute cluster means
    cluster_means = cluster_sums / cluster_counts[:, None]

    # Map back to original points if requested
    if return_mapped:
        mapped_means = cluster_means[inverse_labels]
        if is_1d:
            mapped_means = mapped_means.flatten()
        return mapped_means

    # If original data was 1D, return 1D array
    if is_1d:
        cluster_means = cluster_means.flatten()

    return cluster_means
