
# Module Imports
from . import clustering
from . import manifold
from . import neighbors
from . import plotting
from . import spatial
from . import statistics
from . import utils

# Moduel Aliases
cl = clustering
nn = neighbors
md = manifold
pl = plotting
sp = spatial
stats = statistics

# Surface Function Imports (mostly anndata wrappers)
# from .neighbors._wrapper import neighbors

__all__ = [
    'clustering', 'cl',
    'manifold', 'md',
    'neighbors', 'nn',
    'plotting', 'pl',
    'spatial', 'sp',
    'statistics', 'stats',
    'utils'
]