from scanpy._utils.random import set_igraph_random_state
from igraph import Graph
from typing import List, Dict, Any, Union, Sequence
from numpy import random
import numpy as np
from joblib import Parallel, delayed

scanpy_anyrandom = Union[int, random.RandomState, None]

# ============================================================================ #
#                             Exposed Functions(s)                             #
# ============================================================================ #

def igraph_leiden(
    g: Graph,
    resolution: float = 1.0,
    n_iterations: int = 2,
    random_state: scanpy_anyrandom = 0,
    **clustering_args: Dict[str, Any]
) -> List[int]:
    """
    Run Leiden clustering on a single igraph Graph instance with a fixed random state.

    Parameters
    ----------
    g : Graph
        An igraph Graph object to cluster.
    resolution : float, default=1.0
        Resolution parameter for the Leiden algorithm.
    n_iterations : int, default=2
        Number of iterations to perform in the Leiden algorithm.
    random_state : int, RandomState, or None, default=0
        Random seed or state to ensure reproducible clustering.
    **clustering_args : dict
        Additional keyword arguments to pass to igraph.Graph.community_leiden.

    Returns
    -------
    List[int]
        Cluster membership assignments for each node in the graph.
    """
    clustering_args = dict(clustering_args)
    clustering_args["n_iterations"] = n_iterations
    clustering_args["resolution_parameter"] = resolution
    clustering_args.setdefault("objective_function", "modularity")

    with set_igraph_random_state(random_state):
        part = g.community_leiden(**clustering_args)

    return part.membership

def igraph_leiden_multiseed(
    g: Graph,
    seed_chunk: Sequence[int],
    resolution: float,
    n_iterations: int,
    **kwargs: Dict[str, Any]
) -> List[List[int]]:
    """
    Run igraph_leiden for a list of seeds (chunk) in a single process.

    Parameters
    ----------
    g : Graph
        An igraph Graph object to cluster.
    seed_chunk : Sequence[int]
        A sequence of integer seeds for reproducible clustering.
    resolution : float
        Resolution parameter for the Leiden algorithm.
    n_iterations : int
        Number of iterations per clustering run.
    **kwargs : dict
        Additional keyword arguments passed to igraph_leiden.

    Returns
    -------
    List[List[int]]
        List of cluster membership lists, one per seed.
    """
    return [igraph_leiden(g, resolution, n_iterations, s, **kwargs) for s in seed_chunk]

def igraph_leiden_multiseed_parallel(
    g: Graph,
    resolution: float = 1.0,
    n_iterations: int = 2,
    n_seeds: int = 100,
    n_jobs: int = 8,
    worker_pool: Parallel | None = None,
    **clustering_args: Dict[str, Any]
) -> np.ndarray:
    """
    Run multiple Leiden clustering runs in parallel with different random seeds.

    This function splits the seeds across the requested number of workers, runs
    the clustering in parallel, and returns the results as a NumPy array.

    Parameters
    ----------
    g : Graph
        An igraph Graph object to cluster.
    resolution : float, default=1.0
        Resolution parameter for the Leiden algorithm.
    n_iterations : int, default=2
        Number of iterations per clustering run.
    n_seeds : int, default=100
        Total number of random seeds (i.e., number of clustering runs).
    n_jobs : int, default=8
        Number of parallel workers to use (-1 to use all available CPUs). Only used if worker_pool is None.
    worker_pool: Parallel, Optional
        A pre-initialized joblib multiprocess pool.
    **clustering_args : dict
        Additional keyword arguments passed to igraph_leiden.

    Returns
    -------
    np.ndarray
        Array of shape (n_nodes, n_seeds) containing cluster memberships.
        Each column corresponds to a different seed/run.
    """
    seeds = np.arange(n_seeds)
    seed_chunks = np.array_split(seeds, n_jobs)
    
    worker_pool = worker_pool if worker_pool else Parallel(n_jobs=n_jobs, prefer='processes')
    results = worker_pool(
        delayed(igraph_leiden_multiseed)(
            g=g,
            seed_chunk=seed_chunk,
            resolution=resolution,
            n_iterations=n_iterations,
            **clustering_args
        )
        for seed_chunk in seed_chunks
    )

    # Flatten the list of lists and convert to NumPy array
    return np.array([m for sublist in results for m in sublist], dtype=np.int32).T

# ============================================================================ #
#                             Internal Function(s)                             #
# ============================================================================ #

...

