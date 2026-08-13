from typing import Set, Optional, Callable, Dict, List, Collection
from numpy.typing import NDArray

from collections import defaultdict
from functools import partial
from itertools import combinations

import numpy as np
import json


class ClusterTreeNode :
    def __init__(self, name: str, members: Collection[int], centroid: NDArray[np.float64]) -> None:
        self.name = name
        self.members = set(members)
        self.centroid = centroid

        self.supercluster: Optional[ClusterTreeNode] = None
        self.subclusters: Set['ClusterTreeNode'] = set()
        self.level: int = 0

    @property
    def indexer(self) -> NDArray:
        return np.array(list(self.members))

    def add_child(self, subcluster) -> None:
        assert subcluster.members.issubset(self.members), 'ClusterTree is not hierarchical.'
        subcluster.supercluster = self
        subcluster.level = self.level + 1
        self.subclusters.add(subcluster)
        
    def get_sisters(self) -> list:
        node = self
        sisters = [sister for sister in node.supercluster.subclusters if sister != node]
        while not sisters:
            node = node.supercluster
            sisters = [sister for sister in node.supercluster.subclusters if sister != node]
        
        return sisters
    
    def merge_sister(self, sister) -> None:
        assert sister.level == self.level, 'Only clusters on the same level can be merged.'

        
        # self.name = self.name + '_' + sister.name.split('-')[1]
        self.name = "_".join(sorted(self.name.split('_') + sister.name.split('_')))
        self.members.update(sister.members)
        for subcluster in sister.subclusters:
            subcluster.supercluster = self
            self.subclusters.add(subcluster)

        _n_self = len(self.members)
        _n_other = len(sister.members)
        self.centroid = ((_n_self * self.centroid) + (_n_other * sister.centroid)) / (_n_self + _n_other)

        sister.supercluster.subclusters.remove(sister)
        sister.members.clear()
        sister.subclusters.clear()
        sister.supercluster = None
        
    def get_distance_to(self, other):
        return np.linalg.norm(self.centroid - other.centroid)
    
    def generate_lookup_table(self) -> Dict[str, 'ClusterTreeNode']:
        stack = [self]
        lookup_table = {}

        while stack:
            current = stack.pop()
            lookup_table[current.name] = current
            stack.extend(current.subclusters)

        return lookup_table

    # def generate_lookup_table(self) -> dict:
    #     lookup_table = {self.name: self}

    #     _recurse_lookup_table(self, lookup_table)

    #     return lookup_table

    @classmethod
    def from_labelmat(cls, labelmat: np.ndarray, embedding: np.ndarray, name: Optional[str] = None) -> 'ClusterTreeNode':
        centroid_func = partial(calculate_centroid, embedding=embedding)
        members = np.arange(labelmat.shape[0])
        centroid = centroid_func(members)

        root = cls(name or "root", members, centroid)
        stack = [(root, 0)]

        while stack:
            cluster, depth = stack.pop()
            if depth >= labelmat.shape[1]:
                continue

            level = labelmat[cluster.indexer, depth]
            unique_labels = [label for label in set(level) if isinstance(label, str)]
            for label in unique_labels:
                members = np.where(labelmat[:, depth] == label)[0]

                centroid = centroid_func(members)
                subcluster = ClusterTreeNode(str(label), members, centroid)
                cluster.add_child(subcluster)
                stack.append((subcluster, depth + 1))

        return root

    # @classmethod
    # def from_labelmat(HierTreeCluster, labelmat, embedding, name=None):
    #     centroid_func = partial(calculate_centroid, embedding=embedding)
    #     members = np.arange(labelmat.shape[0])
    #     centroid = centroid_func(members)

    #     # Create the root node
    #     tree = HierTreeCluster(name, members, centroid)

    #     # Build the tree recursively
    #     _recurse_build_tree(labelmat, 0, tree, centroid_func)

    #     return tree
    
    def to_labelmat(self, keep_root: bool = False) -> np.ndarray:
        n_members = len(self.members)
        labelmat = defaultdict(lambda: np.empty(n_members, dtype=object))

        stack = [self]
        while stack:
            cluster = stack.pop()
            labelmat[cluster.level][cluster.indexer] = cluster.name
            stack.extend(cluster.subclusters)

        labelmat = np.vstack([labelmat[level] for level in sorted(labelmat)])
        return labelmat.T if keep_root else labelmat.T[:, 1:]

    # def to_labelmat(self, keep_root=False):
    #     n_members = len(self.members)
    #     labelmat = defaultdict(lambda: np.empty(shape=n_members, dtype=object))

    #     # Recursively populate the labeldict from the tree
    #     _recurse_build_labelmat(self, labelmat)

    #     # Convert the dictionary of levels to a label matrix
    #     labelmat = np.vstack([level.reshape(1, -1) for level in labelmat.values()]).T

    #     if keep_root:
    #         return labelmat
    #     else:
    #         return labelmat[:, 1:]

    @classmethod
    def from_dict(cls, data: Dict) -> 'ClusterTreeNode':
        """
        Creates a ClusterTreeNode object from a dictionary.
        """
        name = data["name"]
        members = data["members"]
        centroid = np.array(data["centroid"], dtype=np.float64)

        cluster = cls(name, members, centroid)
        cluster.level = data["level"]

        # Recursively add subclusters
        for subcluster_data in data["subclusters"]:
            subcluster = cls.from_dict(subcluster_data)
            cluster.add_child(subcluster)

        return cluster

    def to_dict(self) -> Dict:
        """
        Converts the ClusterTreeNode object to a dictionary for JSON/YAML serialization.
        """
        return {
            "name": str(self.name),
            "members": [int(m) for m in self.members],
            "centroid": self.centroid.tolist(),
            "level": int(self.level),
            "subclusters": [subcluster.to_dict() for subcluster in self.subclusters]
        }

    @classmethod
    def from_json(cls, file_path: str) -> 'ClusterTreeNode':
        """
        Loads a ClusterTreeNode object from a JSON file.

        Args:
            file_path (str): The path from which the object will be loaded.

        Returns:
            ClusterTreeNode: The loaded ClusterTreeNode object.
        """
        with open(file_path, 'r') as file:
            data = json.load(file)
            cluster = cls.from_dict(data)
        print(f"ClusterTreeNode '{cluster.name}' loaded from {file_path}")
        return cluster

    def to_json(self, file_path: str) -> None:
        """
        Saves the ClusterTreeNode object to a JSON file.

        Args:
            file_path (str): The path where the object will be saved.
        """
        with open(file_path, 'w') as file:
            json.dump(self.to_dict(), file, indent=4)
        print(f"ClusterTreeNode '{self.name}' saved to {file_path}")

    def to_newick(self, max_depth=None):
        newick = _recurse_build_newick(self, max_depth)
        return f"{newick};"

    def __repr__(self) -> str:
        supercluster = self.supercluster.name if self.supercluster else "None"
        subcluster_names = ', '.join([subcluster.name for subcluster in self.subclusters]) or "None"

        return (f"HierTreeCluster(name='{self.name}', level={self.level}, "
                f"num_members={len(self.members)}, supercluster='{supercluster}', "
                f"subclusters=[{subcluster_names}])")
### Utility Functions
### --------------------------------------------------------------------------------------------------------------------

def calculate_centroid(members, embedding):
    return np.mean(embedding[members, :], axis=0)

### Recursion functions
### --------------------------------------------------------------------------------------------------------------------

# def _recurse_lookup_table(cluster: ClusterTreeNode, lookup_table: dict):
#     for subcluster in cluster.subclusters:
#         lookup_table[subcluster.name] = subcluster
#         _recurse_lookup_table(subcluster, lookup_table)

# def _recurse_build_tree(labelmat: np.ndarray, depth: int, cluster: ClusterTreeNode,
#                        centroid_func: Callable[[set[int]], NDArray[np.float64]]) -> None:
#     level = labelmat[:, depth]
#     # cluster_subset = level[np.array(list(cluster.members))]
#     cluster_subset = level[cluster.indexer]

#     for label in np.unique(cluster_subset):
#         members = np.where(level == label)[0]
#         centroid = centroid_func(members)
#         subcluster = ClusterTreeNode(str(label), set(members), centroid)
#         cluster.add_child(subcluster)

#         if depth + 1 < labelmat.shape[1]:
#             _recurse_build_tree(labelmat, depth + 1, subcluster, centroid_func)

# def _recurse_build_labelmat(cluster: ClusterTreeNode, labelmat: Dict[int, NDArray]):
#     # labelmat[cluster.level][np.array(list(cluster.members))] = cluster.name
#     labelmat[cluster.level][cluster.indexer] = cluster.name

#     for subcluster in cluster.subclusters:
#         _recurse_build_labelmat(subcluster, labelmat)

def _recurse_build_newick(cluster, max_depth=None) -> str:
    # Base case: If the node has no children, it's a leaf node
    if not cluster.subclusters:
        return cluster.name
    
    if max_depth and cluster.level > max_depth:
        return cluster.name

    # Recursive case: Build Newick string for children
    subcluster_newick = ','.join([_recurse_build_newick(subcluster, max_depth) for subcluster in cluster.subclusters])

    return f"({subcluster_newick}){cluster.name}"

### Additional Functions
### --------------------------------------------------------------------------------------------------------------------
def combine_cluster_nodes(clusters: List[ClusterTreeNode], name:str):
    
    assert len(clusters) >= 2, 'Two or more clusters required'
    
    combined_members = set()
    total_members = 0
    weighted_centroid_sum = 0
    
    for cluster in clusters:
        combined_members.update(cluster.members)
        cluster_size = len(cluster.members)
        weighted_centroid_sum += cluster_size * cluster.centroid
        total_members += cluster_size

    # Calculate the new centroid as the weighted average of centroids
    new_centroid = weighted_centroid_sum / total_members

    # Create a new cluster with the combined members and the new centroid
    new_cluster = ClusterTreeNode(name, combined_members, new_centroid)

    return new_cluster

### Aliases
### --------------------------------------------------------------------------------------------------------------------
ClusterTree = ClusterTreeNode