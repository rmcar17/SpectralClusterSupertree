"""Spectral Cluster Supertree.

A scalable and accurate algorithm for merging rooted phylogenetic trees.
"""

from sc_supertree.load import load_trees
from sc_supertree.scs import construct_supertree

__all__ = ["construct_supertree", "load_trees"]
__copyright__ = "Copyright 2023, Robert McArthur"
__license__ = "BSD"
__version__ = "2025.12.9"
