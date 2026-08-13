import numpy as np
from scipy.stats import chi2
from scipy.stats._mannwhitneyu import _broadcast_concatenate
from scipy.stats._stats_py import _rankdata # previously imported from _mannwhitneyu which gets it from here
from typing import Tuple, Optional
from scipy.sparse import spmatrix, issparse

# ============================================================================ #
#                              Exposed Function(s)                             #
# ============================================================================ #

def van_elteren_test(
    X: np.ndarray | spmatrix, 
    Y: np.ndarray | spmatrix, 
    X_strata: np.ndarray, 
    Y_strata: np.ndarray, 
    return_effect_size: bool = False,
    min_stratum_size: int = 10,
    max_stratum_size: int | None = None,
    min_strata: int = 2,
    random_state: int = 0
    ) -> Tuple[np.ndarray, Optional[np.ndarray]]:
    """
    Performs the Van Elteren stratified rank-sum test to compare two groups
    across multiple strata using Mann–Whitney U statistics.

    This implementation supports both dense and sparse matrix inputs and is
    optimized for large-scale biological datasets where sparse expression
    matrices are common.

    Sparse inputs are supported using SciPy sparse matrix formats. In this case,
    data is converted to dense format only after stratum-level subsampling to
    reduce memory overhead and improve computational efficiency.

    This implementation supports:
        - Minimum stratum size filtering for statistical stability
        - Minimum required number of overlapping strata
        - Optional random subsampling of large strata for computational efficiency
        - Sparse matrix inputs (CSR, CSC, COO, and compatible SciPy sparse types).
          CSR highly recommended due to iterative row slicing.

    Parameters
    ----------
    X : np.ndarray or scipy.sparse.spmatrix
        First sample group (n_samples × n_features). May be dense or sparse.

    Y : np.ndarray or scipy.sparse.spmatrix
        Second sample group (n_samples × n_features). May be dense or sparse.

    X_strata : np.ndarray
        Stratum labels corresponding to X samples.

    Y_strata : np.ndarray
        Stratum labels corresponding to Y samples.

    return_effect_size : bool
        If True, also returns stratified effect size estimates.

    min_stratum_size : int
        Minimum number of samples required per group within each stratum.
        Strata with fewer samples in either group are excluded from testing.

    min_strata : int
        Minimum number of valid overlapping strata required for testing.

    max_stratum_size : Optional[int]
        If provided, large strata will be randomly subsampled to this size
        to reduce computational cost and control variance inflation.

        Subsampling is performed prior to dense conversion when sparse matrices
        are provided, ensuring memory efficiency.

    random_state : int
        Random seed controlling subsampling reproducibility.

    Returns
    -------
    pvalues : np.ndarray
        Stratified Van Elteren p-values per feature.

    effect_size : np.ndarray, optional
        Stratified effect size estimates if requested.

    Raises
    ------
    ValueError
        If insufficient overlapping strata are available for testing.

    Notes
    -----
    This implementation approximates the Van Elteren statistic using
    a chi-squared distribution with 1 degree of freedom.

    When max_stratum_size is specified, large strata are randomly
    downsampled using a deterministic pseudorandom generator to ensure
    reproducibility.

    For highly stable inference on very large datasets, consider bootstrap
    aggregation instead of single-pass subsampling.

    Performance is optimized by:
        1. Filtering strata before dense matrix conversion
        2. Subsampling large strata when specified
        3. Vectorized rank aggregation across strata
    """
    assert X.shape[-1] == Y.shape[-1]

    rng = np.random.default_rng(random_state)

    X_is_sparse = issparse(X)
    Y_is_sparse = issparse(Y)
    
    # ------------------------------------------------------------------
    # Phase 1: Validate strata + collect indices
    # ------------------------------------------------------------------
    
    unique_strata = np.unique(X_strata)

    stratum_results = []
    for stratum in unique_strata:
        X_idx = np.where(X_strata == stratum)[0]
        Y_idx = np.where(Y_strata == stratum)[0]
        
        nX = len(X_idx)
        nY = len(Y_idx)
        
        if min(nX, nY) < min_stratum_size:
            continue

        # Subsampling
        if max_stratum_size is not None:
            if len(X_idx) > max_stratum_size:
                X_idx = rng.choice(X_idx, max_stratum_size, replace=False)

            if len(Y_idx) > max_stratum_size:
                Y_idx = rng.choice(Y_idx, max_stratum_size, replace=False)

        stratum_results.append((X_idx, Y_idx))

    # Check overlap power
    if len(stratum_results) < min_strata:
        raise ValueError(
            f"Insufficient overlapping strata. Required {min_strata}, got {len(stratum_results)}"
        )

    # ------------------------------------------------------------------
    # Phase 2 + 3: Compute statistics
    # ------------------------------------------------------------------

    Us = []
    Vs = []
    Ms = []
    Ns = []

    for X_idx, Y_idx in stratum_results:

        sX = X[X_idx].toarray() if X_is_sparse else X[X_idx]
        sY = Y[Y_idx].toarray() if Y_is_sparse else Y[Y_idx]

        U, V, M, N = _van_elteren_core(sX, sY)

        Us.append(U)
        Vs.append(V)
        Ms.append(M)
        Ns.append(N)

    Us = np.array(Us)
    Vs = np.array(Vs)
    Ms = np.array(Ms)[:, None]
    Ns = np.array(Ns)[:, None]

    # ------------------------------------------------------------------
    # Stratified aggregation
    # ------------------------------------------------------------------

    w = (Ms * Ns) / (Ms + Ns + 1e-10)

    den = np.sum(w**2 * Vs, axis=0)
    den[den == 0] = np.nan

    VEs = np.sum(w * Us, axis=0) ** 2 / den
    pvalues = chi2.sf(VEs, df=1)

    # ------------------------------------------------------------------
    # Effect size
    # ------------------------------------------------------------------

    if return_effect_size:

        den = np.sum(w * (Ms + Ns), axis=0)
        den[den == 0] = np.nan

        effect_size = np.sum(w * Us, axis=0) / den

        return effect_size, pvalues

    return pvalues

# ============================================================================ #
#                             Internal Function(s)                             #
# ============================================================================ #

def _van_elteren_core(x: np.ndarray, y: np.ndarray, axis: int = 0) -> Tuple[np.ndarray, np.ndarray, int, int]:
    """
    Performs the core computation for the Van Elteren test, which is based on the Mann-Whitney U statistic,
    and returns the U statistic, variance, and the sizes of the input arrays (M and N).

    The function computes the Mann-Whitney U statistic for two groups (x and y) along the specified axis, 
    taking into account ties in the data, and applies a continuity correction. It also calculates the variance 
    for the statistic, which is corrected for tied ranks.

    Args:
        x (np.ndarray): The first input array (of shape (M, ...)).
        y (np.ndarray): The second input array (of shape (N, ...)).
        axis (int, optional): The axis along which to compute the statistic. Default is 0.

    Returns:
        Tuple[np.ndarray, np.ndarray, int, int]:
            - U: The Mann-Whitney U statistic for the two input arrays.
            - V: The corrected variance of the U statistic.
            - M: The size of the first input array (`x.shape[-1]`).
            - N: The size of the second input array (`y.shape[-1]`).
    """

    # Broadcast and concatenate x and y arrays
    x, y, xy = _broadcast_concatenate(x, y, axis)

    # Get the sizes of the input arrays along the last axis (number of elements)
    M, N = x.shape[-1], y.shape[-1]

    # Rank the concatenated array with ties handled by average ranks
    ranks, t = _rankdata(xy, 'average', return_ties=True)

    # Calculate the Mann-Whitney U statistic for the first group (x)
    R1 = ranks[..., :M].sum(axis=-1)  # Sum of ranks for x
    U = R1 - M * (M + 1) / 2  # Mann-Whitney U statistic, symmetric to U statistic for y

    # Calculate the mean of the U statistic (the expected value of U under the null hypothesis)
    mu = M * N / 2

    # Number of total observations
    n = M + N

    # Calculate the tie term for variance correction
    tie_term = (t**3 - t).sum(axis=-1)

    # Calculate the corrected variance for the U statistic
    V = M * N / 12 * ((n + 1) - tie_term / (n * (n - 1)))  # Tie-corrected variance

    # Center the U statistic around zero (subtract mean U)
    U = U - mu

    # Apply continuity correction (subtract 0.5 depending on the sign of U)
    U -= np.sign(U) * 0.5

    return U, V, M, N