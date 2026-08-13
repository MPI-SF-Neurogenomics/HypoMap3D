import numpy as np

from umap.layouts import optimize_layout_euclidean
from umap.umap_ import make_epochs_per_sample

from hm3d_utils.neighbors import knn_graph, umap_graph
from hm3d_utils.utils import csr_row_normalize
from umap.umap_ import find_ab_params

def project_to_ref_umap(
    ref_X_embed: np.ndarray,
    qry_X_embed: np.ndarray,
    ref_X_umap: np.ndarray,
    spread: float,
    min_dist: float,
    lr: float = 1.0,
    repulsion_strength: float = 1.0,
    negative_sample_rate: float = 5,
    n_epochs: int = 30,
    use_gpu: bool = True,
    verbose: bool = False,
    seed: int = 0,
) -> np.ndarray:
    
    a, b = find_ab_params(spread, min_dist)
    
    # Distances to nearest reference neighbors per query obs
    D_nn = knn_graph(
        reference=ref_X_embed,
        k=15,
        query=qry_X_embed,
        use_gpu=use_gpu
    )

    # Converted to UMAP-style connectivities
    C_nn = umap_graph(
        D_nn,
        k=15
    )
    
    # Create a COO matrix copy
    graph = C_nn.tocoo() # calling graph to match UMAP source code

    # Rowsum softmax from connectivity to probabilty
    csr_row_normalize(C_nn)
    
    # Initial guess of query transformed embedding
    qry_X_umap_init = C_nn @ ref_X_umap
    qry_X_umap_opt = qry_X_umap_init.copy()
    
    # Remove very weak connections
    graph.data[graph.data < (graph.data.max() / float(n_epochs))] = 0.0
    graph.eliminate_zeros()

    # Precompute metrics for optimization
    epochs_per_sample = make_epochs_per_sample(graph.data, n_epochs)
    head = graph.row
    tail = graph.col
    
    qry_X_umap_opt = optimize_layout_euclidean(
        qry_X_umap_opt,
        ref_X_umap.astype(np.float32, copy=True),  # Fixes #179 & #217,
        head,
        tail,
        n_epochs,
        graph.shape[1],
        epochs_per_sample,
        a=a,
        b=b,
        rng_state=np.array([seed, seed, seed], dtype=np.int64),
        gamma=repulsion_strength,
        initial_alpha=lr / 4.0,
        negative_sample_rate=negative_sample_rate,
        parallel=False,
        verbose=verbose,
    )
    
    return qry_X_umap_opt, qry_X_umap_init
