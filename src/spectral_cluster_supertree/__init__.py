"""
Spectral Cluster Supertree

A scalable and accurate algorithm for merging rooted phylogenetic trees.
"""

__all__ = ["load_source_trees", "spectral_cluster_supertree"]
__author__ = "Robert McArthur"
__copyright__ = "Copyright 2024, Robert McArthur"
__credits__ = ["Robert McArthur"]
__license__ = "BSD"
__version__ = "2024.02.26"
__maintainer__ = "Robert McArthur"
__status__ = "Production"

import os

from cogent3 import TreeNode, make_tree

from spectral_cluster_supertree.scs import spectral_cluster_supertree


def load_source_trees(source_tree_file: str | bytes | os.PathLike) -> list[TreeNode]:
    """Load a line-separated file of Newick-formatted trees.

    Parameters
    ----------
    source_tree_file : str | bytes | os.PathLike
        The path to the source tree file.

    Returns
    -------
    list[TreeNode]
        A list of all source trees in the file.
    """
    with open(source_tree_file, "r") as f:
        source_trees: list[TreeNode] = [make_tree(line.strip()) for line in f]  # type: ignore
    return source_trees
