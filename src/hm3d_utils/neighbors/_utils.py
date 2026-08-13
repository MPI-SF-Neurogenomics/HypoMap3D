# This script was discontinued but I am keeping in case I need later

# # Almost entirely ChatGPT written - E

# import numpy as np
# from numba import njit, types, prange

# EPSILON_DBL = 1e-8
# PERPLEXITY_TOLERANCE = 1e-5
# N_STEPS = 100

# @njit(
#     parallel=True,
#     locals={
#         "sum_Pi": types.float64,
#         "sum_disti_Pi": types.float64,
#         "beta": types.float64,
#         "beta_min": types.float64,
#         "beta_max": types.float64,
#         "entropy": types.float64,
#         "entropy_diff": types.float64,
#     },
#     fastmath=True,
# )
# def _knn_bs_perplexity(distances, desired_perplexity):
#     """
#     Compute conditional probabilities for dense kNN distances using squared
#     Euclidean distances for the Gaussian kernel and binary search to match
#     target perplexity.

#     Parameters
#     ----------
#     distances : ndarray of shape (n_samples, n_neighbors)
#         Raw Euclidean distances to neighbors (non-negative).
#     desired_perplexity : float
#         Target perplexity controlling effective number of neighbors.

#     Returns
#     -------
#     P : ndarray of shape (n_samples, n_neighbors)
#         Conditional probability matrix where each row sums to 1. Self-distances
#         (if present on the diagonal) are ignored.
    
#     Notes
#     -----
#     - The binary search is capped at a fixed number of iterations and stops
#       early when the entropy error falls below ``PERPLEXITY_TOLERANCE``.
#     - If the requested perplexity exceeds the number of available neighbors,
#       the entropy constraint may not be exactly achievable.
#     - This function is intended for internal use and mirrors the behavior of
#       perplexity calibration steps in algorithms such as t-SNE.
#     - Zero distances such as self neighbor references are ignored and zeroed.

#     See Also
#     --------
#     _csr_bs_perplexity : Equivalent computation for variable-size neighborhoods
#         stored in CSR sparse matrix format.
#     """
#     n_samples, n_neighbors = distances.shape
#     desired_entropy = np.log(desired_perplexity)
#     P = np.zeros((n_samples, n_neighbors), dtype=np.float64)

#     for i in prange(n_samples):
#         beta_min = -np.inf
#         beta_max = np.inf
#         beta = 1.0

#         for _ in range(N_STEPS):
#             sum_Pi = 0.0
#             # unnormalized probabilities
#             for j in range(n_neighbors):
#                 if distances[i, j] > 0:
#                     dist2 = distances[i, j] ** 2  # square distance for Gaussian
#                     val = np.exp(-dist2 * beta)
#                     P[i, j] = val
#                     sum_Pi += val
#                 else:
#                     P[i, j] = 0.0

#             if sum_Pi == 0.0:
#                 sum_Pi = EPSILON_DBL

#             sum_disti_Pi = 0.0
#             # normalize and compute expected squared distance
#             for j in range(n_neighbors):
#                 pij = P[i, j] / sum_Pi
#                 P[i, j] = pij
#                 sum_disti_Pi += (distances[i, j] ** 2) * pij

#             entropy = np.log(sum_Pi) + beta * sum_disti_Pi
#             entropy_diff = entropy - desired_entropy

#             if abs(entropy_diff) <= PERPLEXITY_TOLERANCE:
#                 break

#             # binary search update
#             if entropy_diff > 0.0:
#                 beta_min = beta
#                 beta = 2.0 * beta if beta_max == np.inf else 0.5 * (beta + beta_max)
#             else:
#                 beta_max = beta
#                 beta = 0.5 * (beta + beta_min) if beta_min != -np.inf else 0.5 * beta

#     return P


# @njit(
#     parallel=True,
#     locals={
#         "sum_Pi": types.float64,
#         "sum_disti_Pi": types.float64,
#         "beta": types.float64,
#         "beta_min": types.float64,
#         "beta_max": types.float64,
#         "entropy": types.float64,
#         "entropy_diff": types.float64,
#         "start": types.int64,
#         "end": types.int64,
#         "nnz_i": types.int64,
#     },
#     fastmath=True,
# )
# def _csr_bs_perplexity(indptr, data, desired_perplexity):
#     """
#     Compute conditional probabilities for sparse CSR distances using squared
#     Euclidean distances for the Gaussian kernel and binary search to match
#     target perplexity.

#     Parameters
#     ----------
#     indptr : ndarray of shape (n_samples + 1,)
#         CSR index pointer array for the sparse distance graph.
#     data : ndarray of shape (nnz,)
#         Raw Euclidean distances of nonzero entries.
#     desired_perplexity : float
#         Target perplexity controlling effective number of neighbors.

#     Returns
#     -------
#     data_P : ndarray of shape (nnz,)
#         Conditional probabilities corresponding to input distances. Each
#         row sums to 1. Zero distances are ignored automatically.
    
#     Notes
#     -----
#     - Each row is processed independently and in parallel using ``prange``.
#     - If a row has fewer nonzero neighbors than the requested perplexity,
#       the entropy constraint may be unattainable.
#     - The function does not modify the sparsity structure; only the values
#       are transformed.
#     - Zero distances such as self neighbor references are ignored and zeroed.
#     - Explicit zeros are generated and must be removed posthoc.

#     See Also
#     --------
#     _knn_bs_perplexity : Dense k-nearest-neighbor version with fixed neighborhood
#         size.
#     """
#     n_samples = indptr.shape[0] - 1
#     desired_entropy = np.log(desired_perplexity)
#     data_P = np.zeros_like(data)

#     for i in prange(n_samples):
#         start = indptr[i]
#         end = indptr[i + 1]
#         nnz_i = end - start
#         if nnz_i == 0:
#             continue

#         beta_min = -np.inf
#         beta_max = np.inf
#         beta = 1.0

#         for _ in range(N_STEPS):
#             sum_Pi = 0.0
#             # unnormalized probabilities
#             for idx in range(start, end):
#                 dist = data[idx]
#                 dist2 = dist ** 2
#                 val = np.exp(-dist2 * beta) if dist > 0 else 0.0
#                 data_P[idx] = val
#                 sum_Pi += val

#             if sum_Pi == 0.0:
#                 sum_Pi = EPSILON_DBL

#             sum_disti_Pi = 0.0
#             # normalize and compute expected squared distance
#             for idx in range(start, end):
#                 pij = data_P[idx] / sum_Pi
#                 data_P[idx] = pij
#                 sum_disti_Pi += (data[idx] ** 2) * pij

#             entropy = np.log(sum_Pi) + beta * sum_disti_Pi
#             entropy_diff = entropy - desired_entropy

#             if abs(entropy_diff) <= PERPLEXITY_TOLERANCE:
#                 break

#             # binary search update
#             if entropy_diff > 0.0:
#                 beta_min = beta
#                 beta = 2.0 * beta if beta_max == np.inf else 0.5 * (beta + beta_max)
#             else:
#                 beta_max = beta
#                 beta = 0.5 * (beta + beta_min) if beta_min != -np.inf else 0.5 * beta

#     return data_P


# # import numpy as np
# # from numba import njit, types, prange

# # EPSILON_DBL = 1e-8
# # PERPLEXITY_TOLERANCE = 1e-5

# # # Specify Numba types for maximum efficiency
# # @njit(
# #     parallel=True,
# #     locals={
# #         "sum_Pi": types.float64,
# #         "sum_disti_Pi": types.float64,
# #         "beta": types.float64,
# #         "beta_min": types.float64,
# #         "beta_max": types.float64,
# #         "entropy": types.float64,
# #         "entropy_diff": types.float64,
# #     },
# #     fastmath=True
# # )
# # def _knn_bs_perplexity(sqdistances, desired_perplexity):
# #     """ [ChatGPT Generated Docstring]
# #     Compute conditional probability distributions with fixed-size neighborhoods
# #     using binary search on the Gaussian bandwidth (precision).

# #     For each sample ``i``, this function finds a precision parameter ``beta_i``
# #     such that the Shannon entropy of the conditional distribution

# #         P(j | i) ∝ exp(-beta_i * sqdistances[i, j])

# #     matches ``log(desired_perplexity)`` up to a specified tolerance. The search
# #     is performed independently for each row using a bounded binary search
# #     procedure.

# #     This implementation assumes a dense k-nearest-neighbor distance matrix
# #     with a uniform number of neighbors per sample and is fully compiled with
# #     Numba for performance.

# #     Parameters
# #     ----------
# #     sqdistances : ndarray of shape (n_samples, n_neighbors)
# #         Squared Euclidean distances to the nearest neighbors for each sample.
# #         Each row corresponds to a sample and contains distances to its
# #         ``n_neighbors`` nearest neighbors. Distances are assumed to be
# #         non-negative.

# #     desired_perplexity : float
# #         Target perplexity for the conditional distributions. The perplexity
# #         controls the effective number of neighbors and is defined as
# #         ``exp(H(P_i))``, where ``H`` is the Shannon entropy.

# #     Returns
# #     -------
# #     P : ndarray of shape (n_samples, n_neighbors)
# #         Conditional probability matrix where each row sums to 1. Entry
# #         ``P[i, j]`` represents the conditional probability of neighbor ``j``
# #         given sample ``i``.

# #     Notes
# #     -----
# #     - The binary search is capped at a fixed number of iterations and stops
# #       early when the entropy error falls below ``PERPLEXITY_TOLERANCE``.
# #     - If the requested perplexity exceeds the number of available neighbors,
# #       the entropy constraint may not be exactly achievable.
# #     - This function is intended for internal use and mirrors the behavior of
# #       perplexity calibration steps in algorithms such as t-SNE.

# #     See Also
# #     --------
# #     _csr_bs_perplexity : Equivalent computation for variable-size neighborhoods
# #         stored in CSR sparse matrix format.
# #     """
# #     n_steps = 100
# #     n_samples, n_neighbors = sqdistances.shape
# #     using_neighbors = n_neighbors < n_samples

# #     desired_entropy = np.log(desired_perplexity)
# #     P = np.zeros((n_samples, n_neighbors), dtype=np.float64)
# #     beta_sum = 0.0

# #     for i in prange(n_samples):
# #         beta_min = -np.inf
# #         beta_max = np.inf
# #         beta = 1.0

# #         for _ in range(n_steps):
# #             sum_Pi = 0.0

# #             # Compute unnormalized probabilities
# #             for j in range(n_neighbors):
# #                 if ((j != i and sqdistances[i, j] > 0) or using_neighbors):
# #                     P[i, j] = np.exp(-sqdistances[i, j] * beta)
# #                     sum_Pi += P[i, j]
# #                 else:
# #                     P[i, j] = 0.0

# #             if sum_Pi == 0.0:
# #                 sum_Pi = EPSILON_DBL

# #             sum_disti_Pi = 0.0
# #             for j in range(n_neighbors):
# #                 P[i, j] /= sum_Pi
# #                 sum_disti_Pi += sqdistances[i, j] * P[i, j]

# #             entropy = np.log(sum_Pi) + beta * sum_disti_Pi
# #             entropy_diff = entropy - desired_entropy

# #             if np.abs(entropy_diff) <= PERPLEXITY_TOLERANCE:
# #                 break

# #             if entropy_diff > 0.0:
# #                 beta_min = beta
# #                 if beta_max == np.inf:
# #                     beta *= 2.0
# #                 else:
# #                     beta = (beta + beta_max) / 2.0
# #             else:
# #                 beta_max = beta
# #                 if beta_min == -np.inf:
# #                     beta /= 2.0
# #                 else:
# #                     beta = (beta + beta_min) / 2.0

# #         beta_sum += beta

# #     return P

# # @njit(
# #     parallel=True,
# #     locals={
# #         "sum_Pi": types.float64,
# #         "sum_disti_Pi": types.float64,
# #         "beta": types.float64,
# #         "beta_min": types.float64,
# #         "beta_max": types.float64,
# #         "entropy": types.float64,
# #         "entropy_diff": types.float64,
# #         "start": types.int64,
# #         "end": types.int64,
# #         "nnz_i": types.int64,
# #     },
# #     fastmath=True,
# # )
# # def _csr_bs_perplexity(indptr, data, desired_perplexity):
# #     """
# #     Compute conditional probability distributions from a CSR sparse distance
# #     matrix using binary search on the Gaussian bandwidth (precision).

# #     This function performs per-row perplexity calibration for a sparse graph
# #     represented in compressed sparse row (CSR) format. For each sample ``i``,
# #     it finds a precision parameter ``beta_i`` such that the Shannon entropy of

# #         P(j | i) ∝ exp(-beta_i * data_ij)

# #     matches ``log(desired_perplexity)`` within a specified tolerance, where
# #     ``data_ij`` are the squared distances to the nonzero neighbors of ``i``.

# #     The sparsity pattern is preserved exactly, and the output probabilities
# #     correspond one-to-one with the input ``data`` array.

# #     Parameters
# #     ----------
# #     indptr : ndarray of shape (n_samples + 1,)
# #         CSR index pointer array. The neighbors of sample ``i`` are stored in
# #         ``data[indptr[i]:indptr[i + 1]]``.

# #     data : ndarray of shape (nnz,)
# #         Squared distances corresponding to the nonzero entries of the sparse
# #         neighborhood graph. Distances are assumed to be non-negative.

# #     desired_perplexity : float
# #         Target perplexity for each conditional distribution. The perplexity
# #         is defined as ``exp(H(P_i))``, where ``H`` is the Shannon entropy.

# #     Returns
# #     -------
# #     data_P : ndarray of shape (nnz,)
# #         Conditional probabilities corresponding to ``data``. For each row
# #         ``i``, the probabilities in ``data_P[indptr[i]:indptr[i + 1]]`` sum
# #         to 1 (up to numerical precision).

# #     Notes
# #     -----
# #     - Each row is processed independently and in parallel using ``prange``.
# #     - If a row has fewer nonzero neighbors than the requested perplexity,
# #       the entropy constraint may be unattainable.
# #     - The function does not modify the sparsity structure; only the values
# #       are transformed.
# #     - This routine is suitable for large, sparse neighborhood graphs as
# #       produced by approximate nearest neighbor search.

# #     See Also
# #     --------
# #     _knn_bs_perplexity : Dense k-nearest-neighbor version with fixed neighborhood
# #         size.
# #     """
# #     n_samples = indptr.shape[0] - 1
# #     desired_entropy = np.log(desired_perplexity)
# #     n_steps = 100

# #     data_P = np.zeros_like(data)

# #     for i in prange(n_samples):

# #         start = indptr[i]
# #         end = indptr[i + 1]
# #         nnz_i = end - start

# #         if nnz_i == 0:
# #             continue

# #         beta_min = -np.inf
# #         beta_max = np.inf
# #         beta = 1.0

# #         for _ in range(n_steps):

# #             sum_Pi = 0.0

# #             # unnormalized probabilities
# #             for idx in range(start, end):
# #                 dist = data[idx]
# #                 if dist > 0.0:
# #                     val = np.exp(-dist * beta)
# #                 else:
# #                     val = 0.0
# #                 data_P[idx] = val
# #                 sum_Pi += val

# #             if sum_Pi == 0.0:
# #                 sum_Pi = EPSILON_DBL

# #             sum_disti_Pi = 0.0

# #             # normalize + expectation
# #             for idx in range(start, end):
# #                 pij = data_P[idx] / sum_Pi
# #                 data_P[idx] = pij
# #                 sum_disti_Pi += data[idx] * pij

# #             entropy = np.log(sum_Pi) + beta * sum_disti_Pi
# #             entropy_diff = entropy - desired_entropy

# #             if abs(entropy_diff) <= PERPLEXITY_TOLERANCE:
# #                 break

# #             if entropy_diff > 0.0:
# #                 beta_min = beta
# #                 if beta_max == np.inf:
# #                     beta *= 2.0
# #                 else:
# #                     beta = 0.5 * (beta + beta_max)
# #             else:
# #                 beta_max = beta
# #                 if beta_min == -np.inf:
# #                     beta *= 0.5
# #                 else:
# #                     beta = 0.5 * (beta + beta_min)

# #     return data_P
