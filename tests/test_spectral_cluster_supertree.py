from pathlib import Path
import pytest
from typing import Sequence
from cogent3 import make_tree, load_tree, TreeNode
from spectral_cluster_supertree import spectral_cluster_supertree

TEST_DATA_DIR = Path("tests/test_data")


def scs_test(in_trees: Sequence[TreeNode], expected: TreeNode):
    result = spectral_cluster_supertree(in_trees).sorted()
    expected = expected.sorted()
    assert result.same_shape(expected)


def load_model_tree_file(model_tree_file):
    return load_tree(TEST_DATA_DIR / model_tree_file)


def load_source_tree_file(source_tree_file):
    with open(TEST_DATA_DIR / source_tree_file) as f:
        source_trees = [make_tree(line.strip()) for line in f]
    return source_trees


def test_agreeable():
    """
    Tests for when spectral clustering is not required.
    """

    tree_1 = make_tree("((a,b),(c,d))")
    tree_2 = make_tree("((a,b),(c,(d,e)))")

    expected = make_tree("((a,b),(c,(d,e)))")

    scs_test([tree_1, tree_2], expected)

    tree_1 = make_tree("(((a,b),(c,d)),(z,(x,y)))")
    tree_2 = make_tree("((a,((f,g),b)),(c,(d,e)))")

    expected = make_tree("(((a,(b,(f,g))),(c,(d,e))),((x,y),z))")

    scs_test([tree_1, tree_2], expected)


@pytest.mark.parametrize(
    "model_tree_file,source_tree_file", [("dcm_model_tree.tre", "dcm_source_trees.tre")]
)
def test_dcm_agreeable(model_tree_file, source_tree_file):
    """
    An example where DCM decomposition means the model tree can always be reproduced.
    """
    model_tree = load_model_tree_file(model_tree_file)
    source_trees = load_source_tree_file(source_tree_file)

    scs_test(source_trees, model_tree)


def test_simple_inconsistency():
    """
    Proper cluster graph shaped:

    a-b-c-d

    Spectral Clustering splits over b/c
    """
    tree_1 = make_tree("(a,(b,c))")
    tree_2 = make_tree("(b,(c,d))")
    tree_3 = make_tree("(d,(a,b))")

    expected = make_tree("((a,b),(c,d))")
    scs_test([tree_1, tree_2, tree_3], expected)


def test_two_squares_inconsitency():
    """
    Proper cluster graph shaped:

    a-b-e-f
    | | | |
    c-d-g-h

    Spectral Clustering splits over b/e d/g

    Subsequent graphs are

    a-b | e-f
        |
    c-d | g-h
    """
    tree_1 = make_tree("((a,b),(c,d))")
    tree_2 = make_tree("((e,f),(g,h))")
    tree_3 = make_tree("(e,(a,c))")
    tree_4 = make_tree("(g,(b,d))")
    tree_5 = make_tree("(a,(e,g))")
    tree_6 = make_tree("(b,(f,h))")
    tree_7 = make_tree("(a,(b,e))")
    tree_8 = make_tree("(h,(d,g))")

    expected = make_tree("(((a,b),(c,d)),((e,f),(g,h)))")
    scs_test([tree_1, tree_2, tree_3, tree_4, tree_5, tree_6, tree_7, tree_8], expected)


def test_simple_contration():
    """
    Small problems where contraction is required
    """
    tree_1 = make_tree("(((a,b),c),(d,e))")
    tree_2 = make_tree("((a,b),(c,d))")

    expected = make_tree("(((a,b),c),(d,e))")
    scs_test([tree_1, tree_2], expected)


def test_size_two_trees():
    """
    Trees with only two leaves. No information gained but a special case nonetheless.
    """
    tree_1 = make_tree("(a,b)")
    tree_2 = make_tree("(b,c)")
    tree_3 = make_tree("(c,d)")
    tree_4 = make_tree("(b,a)")

    expected = make_tree("(a,b,c,d)")
    scs_test([tree_1, tree_2, tree_3], expected)

    scs_test([tree_1], tree_1)
    scs_test([tree_1, tree_1], tree_1)
    scs_test([tree_1, tree_4], tree_1)
