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
