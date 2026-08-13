import numpy as np

from scipy.sparse import csr_matrix

def csr_nz_rowmean(X: csr_matrix) -> np.ndarray:
    "Mean per row of a CSR matrix with respect to only nnz values."
    row_nnz_counts = np.diff(X.indptr)          # number of nnz per row
    row_sums = np.add.reduceat(X.data, X.indptr[:-1])
    row_means = row_sums / row_nnz_counts
    return row_means

def csr_row_normalize(M):
    "Row-wise normalization such that row sums equal 1"
    row_sums = np.array(M.sum(axis=1)).flatten()
    row_sums[row_sums == 0] = 1.0
    
    scale_factors = 1 / row_sums

    M.data *= np.repeat(scale_factors, np.diff(M.indptr))