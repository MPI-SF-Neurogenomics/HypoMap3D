#!/usr/bin/env python3

import numpy as np
import torch

import hm3d_utils

from joblib import Parallel

use_gpu = torch.cuda.is_available()
if use_gpu:
    import rapids_singlecell as rsc
else:
    import scanpy as sc

# ============================================================================ #
#                              Exposed Function(s)                             #
# ============================================================================ #

def run_decomposition(
    adata,
    use_rep: str,
    init_key: str | None = None,
    
    # Decomposition kwargs
    max_levels = 10,
    masking_enabled = True,
    
    # Neighbors kwargs
    n_neighbors: int = 15,
    
    # Leiden kwargs
    resolutions: int | list = [0.1, 0.2, 0.5, 1.0],
    n_seeds: int = 100,
    n_iterations: int = 5,
    n_jobs: int = 8,
    
    # Pruning kwargs
    strata_key: str = 'Sample',
    nn_checks: int = 2,
    min_cluster_size: int = 100,
    max_stratum_size: int = 10_000,
    min_stratum_size: int = 10,
    min_strata: int = 2,
    min_pct_expression: float = 0.1,
    min_specificity: float | int = 0.2,
    significance_threshold: float = 0.001,
    min_marker_genes: int = 20,
    symmetric: bool = True,
    level=0, # set this to initial level (the number of decomp levels that have aready occured)

    # Masking kwargs
    min_nn_agreement: float = 0.5,
    
    random_state=0    
    ) -> np.ndarray:

    resolutions = resolutions if isinstance(resolutions, list) else [resolutions]
    worker_pool = Parallel(n_jobs=n_jobs, prefer='processes')
    
    # If no intial solution provided (not recommended), start fresh from root
    if init_key:
        mask = adata.obs[init_key].isna()
        init_solution = adata.obs[init_key].to_numpy()
        init_solution[mask] = -1
    else:
        init_solution = np.zeros(shape=adata.shape[0])
    labelmat = []
    
    parent_solution = init_solution
    lock = (init_solution == -1) | (init_solution == '-1')
    lock_active = False
    while not np.all(lock):
        
        labelmat.append(parent_solution)
        
        if level >= max_levels:
            break
        
        solution = np.full_like(parent_solution, fill_value=-1)
        
        for i, cluster in enumerate(set(parent_solution[~lock])):
            no_decomposition = False
            
            parent_mask = parent_solution == cluster
            
            subset = adata[parent_mask]
            
            # 1. Compute Consensus Clusters ---------------------------------------------- #
            if use_gpu:
                rsc.pp.neighbors(
                    subset, 
                    n_neighbors=n_neighbors,
                    use_rep=use_rep
                )
            else:
                sc.pp.neighbors(
                    subset,
                    n_neighbors=n_neighbors,
                    use_rep=use_rep
                )


            connectivities_key = subset.uns['neighbors']['connectivities_key']
            G = hm3d_utils.utils.get_igraph_from_adjacency(subset.obsp[connectivities_key])
            
            if level < len(resolutions):
                r = resolutions[level]
                # keep lock disabled as higher resolutions may result in decomposition, where lower resolutions didn't
            else:
                r = resolutions[-1]
                lock_active=True # fixed resolutions, enable locking
            
            solutions = hm3d_utils.cl.igraph_leiden_multiseed_parallel(
                g=G,
                resolution=r,
                n_iterations=n_iterations,
                n_seeds=n_seeds,
                worker_pool=worker_pool
            )
            
            nunique = np.array([len(set(solution)) for solution in solutions.T])
            single_mask = nunique == 1
            
            # If majority of solutions don't decompose, then path terminates and skip
            fraction_single = np.mean(single_mask)
            if fraction_single < 0.5:
                solutions = solutions[:, ~single_mask]
                try:
                    consensus, _ = hm3d_utils.cl.find_consensus(
                        solutions=solutions,
                        random_state=random_state
                    )
                except np.linalg.LinAlgError:
                    no_decomposition = True
            else:
                no_decomposition = True
            
            
            # 2. Prune Consensus Clusters ------------------------------------------------ #
            
            if not no_decomposition:
             
                embedding = subset.obsm[use_rep]
                
                pruned = hm3d_utils.cl.prune_clusters(
                    membership=consensus,
                    gexpr=subset.X,
                    embedding=embedding,
                    strata=subset.obs[strata_key].values,
                    nn_checks=nn_checks,
                    min_cluster_size=min_cluster_size,
                    max_stratum_size=max_stratum_size,
                    min_stratum_size=min_stratum_size,
                    min_strata=min_strata,
                    min_pct_expression=min_pct_expression,
                    min_specificity=min_specificity,
                    significance_threshold=significance_threshold,
                    min_marker_genes=min_marker_genes,
                    symmetric=symmetric,
                    random_state=random_state,
                )
                
                # Unique ids for each new cluster across all levels
                pruned = _cantor_pairing(level, _cantor_pairing(i, pruned))
                
                if len(set(pruned)) == 1:
                    no_decomposition = True
                
                elif masking_enabled:
                    # 3. Mask noisy elements (optional) ------------------------------------------ #
                    # To stop noisy error in element assignment propagating we drop them for later
                    
                    distances_key = subset.uns['neighbors']['distances_key']
                    D = subset.obsp[distances_key]
                    
                    nn_agreement = hm3d_utils.stats.neighborhood_agreement(
                        D=D,
                        memberships=pruned
                    )
                    uncertain_mask = nn_agreement < min_nn_agreement
                    pruned[uncertain_mask] = -1
            
            # 4. Propagate or lock when no decomposition occurs -------------------------- #
            if no_decomposition:
                if lock_active:
                    lock[parent_mask] = True # lock branch, propagate downstream
                else:
                    solution[parent_mask] = parent_solution[parent_mask] # propagate immediately
                continue            
            else:
                solution[parent_mask] = pruned
        
        # 5. Propagate unmerged via lock when active --------------------------------- #
        if lock_active:
            solution[lock] = parent_solution[lock]
        
        # Update lock to True for masked values (independant of lock activity)
        lock[solution == -1] = True
        
        # 6. Write solution to labelmat ---------------------------------------------- #
        parent_solution = solution
        
        print(f'Step finished with {np.mean(lock) * 100: 0.2f}% of obs locked and {len(set(parent_solution[parent_solution!=-1]))} clusters.', ' '*90)
        
        level +=1

    labelmat = np.asarray(labelmat, dtype=str).T

    return labelmat

# ============================================================================ #
#                               Helper Functions                               #
# ============================================================================ #

def _cantor_pairing(a, b_arr):
    a = np.int64(a)
    b_arr = np.asarray(b_arr, dtype=np.int64)
    s = a + b_arr
    return (s * (s + 1)) // 2 + b_arr