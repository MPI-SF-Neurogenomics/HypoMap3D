import numpy as np
from numba import njit

# ============================================================================ #
#                              Exposed Function(s)                             #
# ============================================================================ #

def specificity(
    logFC, 
    pctA,
    pctB,
    keep_sign=False
):
    pass

def log2_foldchange(
    meanA,
    meanB,
    log1p=False   
):
    """_summary_

    Parameters
    ----------
    meanA : _type_
        _description_
    meanB : _type_
        _description_
    log1p : bool, optional
        _description_, by default False
    
    Notes
    -----
    - Seurat style logFC computation via pseudocount 1/n_cells
    """
    pass

# @njit
# def pct_mean_gexpr(
#     data: np.ndarray[float], 
#     indices: np.ndarray[int], 
#     n_obs: int, 
#     n_genes: int, 
#     threshold: int | float = 0
# ) -> tuple[np.ndarray, np.ndarray]:
#     """
#     Compute per-gene mean expression and fraction of cells expressing above a threshold
#     from a sparse matrix represented in CSR-like format (nonzero values and column indices).

#     This function performs a single-pass, branchless accumulation over nonzero entries
#     for maximal performance with Numba JIT compilation.

#     Parameters
#     ----------
#     data : np.ndarray[float]
#         1D array of nonzero values from a sparse matrix (typically `adata.X.data`).
#     indices : np.ndarray[int]
#         Column indices corresponding to each entry in `data` (typically `adata.X.indices`).
#     n_obs : int
#         Total number of observations (rows, e.g., cells) in the original matrix.
#     n_genes : int
#         Total number of features (columns, e.g., genes) in the original matrix.
#     threshold : int or float, default=0
#         Value above which a gene is considered “expressed” in a given observation.
#         Must be >= 0 for sparse matrices where zeros are implicit.

#     Returns
#     -------
#     pct_gexpr : np.ndarray[float]
#         Fraction of observations where each gene is expressed above `threshold`.
#         Shape: `(n_genes,)`.
#     mean_gexpr : np.ndarray[float]
#         Mean expression of each gene across all observations, including zeros.
#         Shape: `(n_genes,)`.

#     Notes
#     -----
#     - Operates efficiently on CSR-like sparse matrices without densifying.
#     """
    
#     # allocate outputs
#     _colsums = np.zeros(n_genes, dtype=np.float32)
#     _counts  = np.zeros(n_genes, dtype=np.int32)

#     # local aliases for fast access
#     _cs  = _colsums
#     _ct  = _counts
#     _dat = data
#     _ind = indices
#     _thresh = threshold

#     _nnz = len(_dat)

#     # single linear pass
#     for _i in range(_nnz):
#         _v = _dat[_i]
#         _j = _ind[_i]
#         _cs[_j] += _v
#         _ct[_j] += _v > _thresh  # branchless

#     _inv_n_obs = 1.0 / n_obs
#     mean_gexpr = _cs * _inv_n_obs
#     pct_gexpr  = _ct * _inv_n_obs

#     return pct_gexpr, mean_gexpr