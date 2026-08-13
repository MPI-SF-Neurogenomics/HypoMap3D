import warnings

import anndata as ad
import numpy as np

from scipy.sparse import issparse

def view_to_actual(adata: ad.AnnData) -> ad.AnnData:
    if adata.is_view:
        adata._init_as_actual(adata.copy())
    
    return adata

def feature_percentages(
    adata: ad.AnnData,
    min_value: int | float = 0
    ) -> np.ndarray:
    """
    Compute the fraction of cells expressing each feature above a threshold.

    This function calculates, for each feature (gene, protein, or other
    measured variable) in an AnnData object, the proportion of cells
    where the feature's value exceeds a given minimum threshold.

    Parameters
    ----------
    adata : anndata.AnnData
        The AnnData object containing the feature matrix `.X`. It can be
        dense or sparse (CSR) format. Features are assumed to be stored
        in columns.
    min_value : float, default=0
        The minimum value for a feature to be considered "expressed"
        or "present" in a cell. Values strictly greater than `min_value`
        are counted.

    Returns
    -------
    numpy.ndarray
        A 1D array of length equal to the number of features, containing
        the fraction of cells in which each feature exceeds `min_value`.

    Notes
    -----
    - This function is optimized for sparse matrices and avoids creating
      dense intermediate arrays.
    - For a threshold of 0, it returns the standard "detection rate"
      for each feature.
    - The returned values are in the range [0, 1].
    """
    n_cells, n_features = adata.shape
    mask = adata.X.data > min_value
    counts = np.bincount(adata.X.indices[mask], minlength=n_features)
    return counts / n_cells

def write_clean_h5ad(
    adata: ad.AnnData, 
    outfile: str, 
    **kwargs
    ) -> None:
    """
    Clean and write an AnnData object to disk in .h5ad format with memory and consistency optimizations.

    This function performs lightweight preprocessing on the AnnData object prior to saving,
    primarily aimed at reducing file size and avoiding common serialization issues in large
    single-cell datasets.

    Steps performed:
    - Downcasts `adata.X` from float64 to float32 (if applicable) to reduce memory footprint.
    - Removes explicit zero entries from sparse matrices.
    - Downcasts sparse index arrays to int32 when safe (for storage efficiency).
    - Removes color-related metadata entries from `adata.uns`.
    - Issues a warning if `'log1p'` is found in `.uns`, which may indicate normalized data.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix to be cleaned and written.
    outfile : str or Path
        Destination path for the output `.h5ad` file.
    **kwargs
        Additional keyword arguments passed to `adata.write_h5ad()`.

    Returns
    -------
    None
        The function writes the file to disk and does not return a value.

    Notes
    -----
    This function modifies the input AnnData object in place. If preservation of the original
    object is required, pass a copy (e.g., `adata.copy()`).

    The optimizations performed here prioritize storage efficiency and compatibility, but may
    remove some metadata (e.g., color annotations) that are useful for visualization.
    
    The downcast of `adata.X.indices` to `np.int32` is reverted whenever `.X` is copied. This patch
    can also result in errors such as calling `adata.X.eliminate_zeros`.

    Warnings
    --------
    A warning is raised if `'log1p'` is detected in `.uns`, as this may indicate that the data
    stored in `.X` has already been log-transformed or normalized.

    """
    
    # Downcast float64 to float32 for size reduction
    if adata.X.dtype == np.float64:
        adata.X = adata.X.astype(np.float32)
    
    # Eliminate potential zeros
    if issparse(adata.X):
        try:
            adata.X.eliminate_zeros()
        except:
            warnings.warn("Zeros already eliminated")
    
    # dtype conversion trick to reduce size
    if np.max(adata.X.indices) < 2**31:
        adata.X.indices = adata.X.indices.astype(np.int32, copy=False)
    
    # remove added color information to keys and warn if counts are normalized
    for key in adata.uns.keys():
        if '_color' in key:
            del adata.uns[key]
        
        if key == 'log1p':
            warnings.warn('`log1p` found in adata.uns. Are you sure you want to save normalized counts?')
    
    adata.write_h5ad(outfile, **kwargs)