from pathlib import Path


import pytest
from cogent3 import make_tree
from spectral_cluster_supertree import spectral_cluster_supertree
from spectral_cluster_supertree.scs.scs import tuple_sorted


def scs_test(in_trees, expected):
    result = spectral_cluster_supertree(in_trees).sorted()
    assert str(result) == str(expected)


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

def test_tuple_sorted():
    assert tuple_sorted([3,1,2]) == (1,2,3)
