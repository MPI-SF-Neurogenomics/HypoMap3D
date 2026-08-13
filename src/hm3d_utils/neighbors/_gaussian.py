import numpy as np
from numba import njit, types, prange
from scipy.sparse import csr_matrix

from typing import Literal, Callable

# ============================================================================ #
#                              Exposed Function(s)                             #
# ============================================================================ #

def gauss_graph(
    D: csr_matrix,
    k: int | float,
    sq_distances: bool = True,
    sigma: Literal['median', 'std'] = 'std'
) -> csr_matrix:
    ...


def gauss_knn(
    knn_dists: np.ndarray,
    k: int | float,
    sq_distances: bool = True,
    sigma: Literal['median', 'std'] = 'std'
) -> np.ndarray:
    ...

# ============================================================================ #
#                              Helper Function(s)                              #
# ============================================================================ #

@njit(parallel=True)
def _gauss_csr(
    distances: np.ndarray,
    indptr: np.ndarray,
    sig_func: Callable
) -> np.ndarray:
    ...

@njit(parallel=True)
def _gauss_knn(
    knn_dists: np.ndarray,
    sig_func: Callable
) -> np.ndarray:
    ...


@njit(fastmath=True)
def _gauss_kernel_transform(
    distances: np.ndarray,
    sigma: int | float
) -> np.ndarray:
    ...