import os
from pathlib import Path

from cogent3 import PhyloNode, make_tree


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
