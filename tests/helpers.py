from pathlib import Path
from typing import Literal, Sequence
from cogent3 import TreeNode, load_tree, make_tree

from spectral_cluster_supertree.scs import spectral_cluster_supertree

TEST_DATA_DIR = Path("tests/test_data")


def scs_test(
    in_trees: Sequence[TreeNode],
    expected: TreeNode,
    weights: Sequence[float] | None = None,
    pcg_weighting: Literal["one", "branch", "depth"] = "one",
    contract_edges: bool = True,
):
    result = spectral_cluster_supertree(
        in_trees,
        weights=weights,
        pcg_weighting=pcg_weighting,
        contract_edges=contract_edges,
    ).sorted()
    expected = expected.sorted()
    assert result.same_shape(expected), str(result) + " != " + str(expected)


def load_expected_tree_file(model_tree_file):
    return load_tree(TEST_DATA_DIR / model_tree_file)


def load_source_tree_file(source_tree_file):
    with open(TEST_DATA_DIR / source_tree_file) as f:
        source_trees = [make_tree(line.strip()) for line in f]
    return source_trees
