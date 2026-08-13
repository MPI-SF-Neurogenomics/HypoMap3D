import numpy as np
import pandas as pd

from ..utils._anndata import feature_percentages
from ._van_elteren import van_elteren_test

def find_markers(
    adata_A,
    adata_B,
    strata_key,
    min_pct_expression=0.1,
    min_specificity=0.3,
    min_stratum_size=10,
    max_stratum_size=10_000,
    min_strata=2,
    random_state=0,
):
    
    n_A, _ = adata_A.shape
    n_B, _ = adata_B.shape

    pct_A = feature_percentages(adata_A, min_value=0)
    pct_B = feature_percentages(adata_B, min_value=0)
    pct_mask = (pct_A >= min_pct_expression) | (pct_B >= min_pct_expression)

    # Assumed expression is not log1p normalized
    mean_A = adata_A.X.mean(axis=0).A1
    mean_B = adata_B.X.mean(axis=0).A1
    log2_foldchange = np.log2(mean_A + (1/n_A)) - np.log2(mean_B + (1/n_B)) # seurat-style
    marker_sign_vec = np.sign(log2_foldchange)

    rel_pct_diff = np.abs(pct_A - pct_B) / (pct_A + pct_B + 1e-10)
    specificity = rel_pct_diff * np.abs(log2_foldchange)

    specificity_mask = specificity >= min_specificity

    # Gene to be tested
    valid_genes = pct_mask & specificity_mask

    X_A = adata_A.X[:, valid_genes]
    X_B = adata_B.X[:, valid_genes]
    
    labels_A = adata_A.obs[strata_key].to_numpy()
    labels_B = adata_B.obs[strata_key].to_numpy()

    _, inverse = np.unique(np.concatenate([labels_A, labels_B]), return_inverse=True)

    strata_A = inverse[:len(labels_A)]
    strata_B = inverse[len(labels_A):]

    effect_sizes, pvalues = van_elteren_test(
        X=X_A, 
        Y=X_B, 
        X_strata=strata_A, 
        Y_strata=strata_B,
        min_stratum_size=min_stratum_size,
        max_stratum_size=max_stratum_size,
        min_strata=min_strata,
        return_effect_size=True,
        random_state=random_state
    )

    padj = pvalues * np.sum(valid_genes)

    direction_map = {1: 'A', 0: pd.NA, -1: 'B'}

    results = pd.DataFrame(
        {
            'direction': [direction_map[x] for x in marker_sign_vec[valid_genes]],
            'specificity': specificity[valid_genes],
            'pctA': pct_A[valid_genes],
            'pctB': pct_B[valid_genes],
            'meanA': mean_A[valid_genes],
            'meanB': mean_B[valid_genes],
            'log2FC': log2_foldchange[valid_genes],
            'p': pvalues,
            'padj': np.clip(padj, 0, 1),
            'VE': effect_sizes
        },
        index = adata_A.var_names[valid_genes]
    )
    
    results['rank'] = results['specificity'].copy()
    results.loc[results['direction'] == 'A', 'rank'] = results['pctA'] * results['specificity']
    results.loc[results['direction'] == 'B', 'rank'] = results['pctB'] * results['specificity']
    results['rank'] = results['rank'].rank(method='min', ascending=False).astype(int)
    
    results = results.sort_values('rank')

    return results