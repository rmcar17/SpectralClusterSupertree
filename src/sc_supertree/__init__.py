"""Spectral Cluster Supertree.

A scalable and accurate algorithm for merging rooted phylogenetic trees.
"""

__all__ = ["construct_supertree", "load_trees"]
__copyright__ = "Copyright 2023, Robert McArthur"
__license__ = "BSD"
__version__ = "2025.6.11"

import os
from pathlib import Path

from cogent3 import PhyloNode, make_tree

from sc_supertree.scs import construct_supertree


def load_trees(source_tree_file: str | os.PathLike) -> list[PhyloNode]:
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
    with Path(source_tree_file).open() as f:
        source_trees: list[PhyloNode] = [make_tree(line.strip()) for line in f]
    return source_trees
