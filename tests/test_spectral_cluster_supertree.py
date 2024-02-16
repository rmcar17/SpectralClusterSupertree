from cogent3 import make_tree
from spectral_cluster_supertree import spectral_cluster_supertree


def scs_test(in_trees, expected):
    result = spectral_cluster_supertree(in_trees).sorted()
    assert result.same_shape(expected)


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


def test_simple_inconsistency():
    tree_1 = make_tree("(a,(b,c))")
    tree_2 = make_tree("(b,(c,d))")
    tree_3 = make_tree("(d,(a,b))")

    expected = make_tree("((a,b),(c,d))")
    scs_test([tree_1, tree_2, tree_3], expected)
