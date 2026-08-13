import numpy as np

from scipy.sparse import diags, linalg
from sklearn.cluster._spectral import discretize
from sklearn.preprocessing import OneHotEncoder

from ..statistics import element_agreement

# ---------------------------------------------------------------------------- #
#                              Exposed Function(s)                             #
# ---------------------------------------------------------------------------- #

def find_consensus(
    solutions: np.ndarray,
    random_state: int = 0,
) -> tuple[np.ndarray, np.ndarray]:
    """ [ChatGPT + Manual Review]
    Compute a consensus clustering from multiple clustering solutions using Hybrid Bipartite Graph
    Formulation (HBGF).

    This function takes multiple clustering solutions for the same set of elements (e.g., cells),
    constructs a sparse adjacency matrix, computes a normalized graph Laplacian, and then performs
    truncated singular value decomposition (SVD). A discrete partition is obtained using the
    spectral embedding, and the consensus labels are chosen to maximize agreement across all input
    solutions using elem_sim_scores (`hm3d_utils.stats.element_agreement`).

    Parameters
    ----------
    solutions : np.ndarray, shape (n_elements, n_solutions)
        Array containing multiple clustering solutions. Each row corresponds to an element (e.g., cell),
        and each column corresponds to a clustering solution. Each entry is an integer cluster label.
    random_state : int, default=0
        Seed for reproducibility of the SVD and discretization steps.

    Returns
    -------
    best_labels : np.ndarray, shape (n_elements,)
        Array of consensus cluster labels for each element. These labels maximize agreement across
        the input clustering solutions.
    best_agreement : np.ndarray, shape (n_solutions,)
        Array of mean agreement scores between the consensus labels and each input clustering solution.
        Each entry corresponds to the agreement of the consensus all input solutions.

    Notes
    -----
    - The function builds a one-hot-like sparse adjacency matrix from the input solutions,
      where columns represent all clusters across all solutions.
    - A normalized graph Laplacian is computed as `L = A * D_inv` where `D_inv` is the inverse
      square root of the column sums of `A`.
    - SVD is performed on `L` using either 'arpack' or 'propack', depending on matrix size.
    - Left and right singular vectors are stacked and then discretized to produce candidate labels.
    - The consensus labels are selected as those maximizing the mean element-wise agreement
      across all input solutions.
    - Any NaNs in `solutions` should be handled prior to calling this function, as they can affect
      cluster counting and adjacency encoding.
    """  
    
    n_clusters_vec = np.apply_along_axis(lambda sol: len(np.unique(sol)), arr=solutions, axis=0)
    k_max = np.max(n_clusters_vec)
    
    # A = _construct_adj_matrix(solutions)
    ohe = OneHotEncoder(dtype=bool)
    A = ohe.fit_transform(solutions)
    A = A.tocsc()

    # Compute graph laplacian matrix
    degree = np.array(A.sum(axis=0)).ravel()
    degree[degree == 0] = 1 # avoid potential zero division (v. unlikely)
    D_inv = diags(1 / np.sqrt(degree))
    L = A @ D_inv
    
    # row-oriented matrix -> speed ups SVD to convert
    L = L.tocsr()

    # Compute left and right singular values
    svd_k = min(min(L.shape)-1, k_max)
    solver = 'propack' if (A.shape[0] > 350_000 or A.shape[1] > 50_000) else 'arpack'
    U, _, V = linalg.svds(L, k=svd_k, solver=solver, which='LM', random_state=random_state)

    U = U / np.linalg.norm(U, axis=1, keepdims=True)
    V = V / np.linalg.norm(V, axis=0, keepdims=True)

    n = U.shape[0]
    UVt = np.vstack([U, V.T])

    best = -1
    best_labels = None

    # Find k clusters that maximize stability
    for k in np.unique(n_clusters_vec):

        UVt_sliced = UVt[:, :k]

        labels = discretize(UVt_sliced, random_state=random_state)
        labels = labels[:n] # clip cluster node memberships

        agreement = element_agreement(labels, solutions.T)
        score = np.mean(agreement)

        if score > best:
            best = score
            best_labels = labels
            best_agreement = agreement
    
    return best_labels, best_agreement

# ---------------------------------------------------------------------------- #
#                            Internal Function(s)                              #
# ---------------------------------------------------------------------------- #

# def _construct_adj_matrix(
#     matrix
# ):
#     """
#     Construct a sparse adjacency matrix from a clustering matrix.

#     This function constructs a sparse adjacency matrix from a clustering matrix, where each column represents
#     a clustering solution, and each row represents an element being clustered. It converts each clustering
#     solution into a binary matrix and concatenates them horizontally to form the adjacency matrix.

#     Parameters
#     ----------
#     matrix : ndarray
#         A 2D array where each column represents a clustering solution, and each row represents an element being clustered.

#     Returns
#     -------
#     csc_matrix
#         A sparse adjacency matrix in Compressed Sparse Column (CSC) format.
#     """
    
#     n, m = matrix.shape
#     encoded = np.empty_like(matrix, dtype=np.int64)

#     offset = 0
#     for i in range(m):
#         solution = matrix[:, i]
#         unique_labels, encoded_col = np.unique(solution, return_inverse=True)
#         encoded[:, i] = encoded_col + offset
#         offset += len(unique_labels)
        
#     A_rows, A_cols, A_data = _compute_adjacency_coords(encoded)

#     A = csc_matrix((A_data, (A_rows, A_cols)),
#                    shape=(n, offset),
#                    dtype=np.bool_)
    
#     return A

# @njit
# def _compute_adjacency_coords(matrix):
#     n, m = matrix.shape
#     nnz = n * m

#     rows = np.empty(nnz, dtype=np.int64)
#     cols = np.empty(nnz, dtype=np.int64)
#     data = np.ones(nnz, dtype=np.bool_)

#     idx = 0
#     for j in range(m):
#         sol = matrix[:, j]
#         for i in range(n):
#             rows[idx] = i
#             cols[idx] = sol[i]
#             idx += 1

#     return rows, cols, data