from typing import Iterable, Sequence, Optional, Set
from cogent3.core.tree import TreeNode, TreeBuilder,PhyloNode


def spectral_cluster_supertree(
    trees: Sequence[TreeNode], weights: Optional[Sequence[float]]
) -> TreeNode:
    """_summary_

    Args:
        trees (Iterable[TreeNode]): _description_
        weights (Optional[float]): _description_

    Returns:
        TreeNode: _description_
    """

    assert len(trees) >= 1, "there must be at least one tree"

    if weights is None:
        weights = [1.0 for _ in range(len(trees))]

    assert len(trees) == len(weights), "trees and weights must be of same length"

    all_names = get_all_tip_names(trees)

    if len(all_names) <= 2:


    return TreeNode()


def get_all_tip_names(trees: Iterable[TreeNode]) -> Set:
    names = set()
    for tree in trees:
        names.update(tree.get_tip_names())
    return names

def tip_names_to_tree(tip_names: Iterable) -> TreeNode:
    tree_builder = TreeBuilder(constructor=TreeNode).create_edge # Incorrectly causes "type error" on input TreeNode due to type detection system
    tips = [tree_builder([], tip_name, {}) for tip_name in tip_names]
    tree = tree_builder(tips, "root", {}) # Might need to change "root" to something else, it is likely only temporarily the root.
    return tree



