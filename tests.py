from cogent3 import make_tree
from spectral_cluster_supertree import spectral_cluster_supertree


# print(spectral_cluster_supertree([make_tree("(a)")]))
def scs_test(in_trees, expected, verbose=False):
    result = spectral_cluster_supertree(in_trees)
    if verbose:
        print("In:", in_trees)
        print("Out:", result)
    print(result, expected)
    # TODO: Should this be result.same_topology?
    # assert result == expected


def test_agreeable():
    """
    Tests for when spectral clustering is not required.
    """

    tree1 = make_tree("((a,b),(c,d))")
    tree2 = make_tree("((a,b),(c,(d,e)))")

    expected = make_tree("((a,b),(c,(d,e)))")

    scs_test([tree1, tree2], expected, verbose=True)

    tree1 = make_tree("(((a,b),(c,d)),(z,(x,y)))")
    tree2 = make_tree("((a,((f,g),b)),(c,(d,e)))")

    expected = make_tree("(((a,(b,(f,g))),(c,(d,e))),((x,y),z))")

    scs_test([tree1, tree2], expected, verbose=True)


if __name__ == "__main__":
    test_agreeable()

    a = make_tree("(((a,b),c),(d,e))")
    b = make_tree("((a,b),(c,d))")
    c = make_tree("((a,e),(c,d))")
    d = make_tree("((b,e),(c,d))")
    print(
        spectral_cluster_supertree([a, b], [1, 1])
    )  # Interesting results for this test case, around 6 for the second
    # print(spectral_cluster_supertree([a, b, c]))
    # print(spectral_cluster_supertree([a, b, c, d]))
