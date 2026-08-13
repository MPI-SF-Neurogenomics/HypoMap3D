import igraph as ig
import numpy as np

from scipy.sparse import spmatrix

def get_igraph_from_adjacency(adjacency: spmatrix | np.ndarray, *, directed: bool = False) -> ig.Graph:
    """Convert adjacency matrix to igraph.Graph."""
    sources, targets = adjacency.nonzero()
    g = ig.Graph(directed=directed)
    g.add_vertices(adjacency.shape[0])
    g.add_edges(zip(sources, targets))
    
    # add weights if present
    if np.issubdtype(adjacency.dtype, np.number):
        g.es["weight"] = adjacency[sources, targets].A1  # .A1 = flatten
    
    return g