from collections.abc import Sequence
from typing import Any, Literal

import pytest
from cogent3 import PhyloNode, make_tree
from cogent3.app.composable import NotCompleted
from helpers import load_expected_tree_file, load_source_tree_file

from sc_supertree import construct_supertree


def scs_test(
    in_trees: Sequence[PhyloNode],
    expected: PhyloNode,
    weights: Sequence[float] | None = None,
    pcg_weighting: Literal["one", "branch", "depth", "bootstrap"] = "one",
    *,
    contract_edges: bool = True,
) -> None:
    result = construct_supertree(
        in_trees,
        weights=weights,
        pcg_weighting=pcg_weighting,
        contract_edges=contract_edges,
    ).sorted()
    expected = expected.sorted()
    assert result.same_shape(expected), f"{result} != {expected}"


def test_agreeable() -> None:
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
    ("model_tree_file", "source_tree_file"),
    [("dcm_model_tree.tre", "dcm_source_trees.tre")],
)
def test_dcm_agreeable(model_tree_file: str, source_tree_file: str) -> None:
    """
    An example where DCM decomposition means the model tree can always be reproduced.
    """
    model_tree = load_expected_tree_file(model_tree_file)
    source_trees = load_source_tree_file(source_tree_file)

    scs_test(source_trees, model_tree)


def test_simple_inconsistency() -> None:
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


def test_two_squares_inconsitency() -> None:
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


def test_simple_contration() -> None:
    """
    Small problems where contraction is required
    """
    tree_1 = make_tree("(((a,b),c),(d,e))")
    tree_2 = make_tree("((a,b),(c,d))")

    expected = make_tree("(((a,b),c),(d,e))")
    scs_test([tree_1, tree_2], expected)


def test_size_two_trees() -> None:
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


def test_simple_weights() -> None:
    """
    Tests tree weighting works as expected
    """
    tree_1 = make_tree("(a,(b,c))")
    tree_2 = make_tree("(c,(a,b))")

    scs_test([tree_1, tree_2], tree_1, weights=[2, 1])
    scs_test([tree_1, tree_2], tree_1, weights=[1.001, 1])

    scs_test([tree_1, tree_2], tree_2, weights=[1, 2])
    scs_test([tree_1, tree_2], tree_2, weights=[1, 1.001])


def test_depth_weighting() -> None:
    """
    Test that depth pcg weighting works appropriately
    (branch equivalent when no lengths).
    """
    tree_1 = make_tree("(a,(b,(c,(d,e))))")
    tree_2 = make_tree("(d,(f,(a,b)))")

    expected_one = make_tree("((f,a),(b,(c,(d,e))))")
    expected_depth = make_tree("((f,(a,b)),(c,(d,e)))")

    scs_test(
        [tree_1, tree_2],
        expected_one,
        pcg_weighting="one",
        contract_edges=False,
    )
    scs_test(
        [tree_1, tree_2],
        expected_depth,
        pcg_weighting="depth",
        contract_edges=False,
    )
    scs_test(
        [tree_1, tree_2],
        expected_depth,
        pcg_weighting="branch",
        contract_edges=False,
    )


def test_branch_weighting() -> None:
    """
    Test branch lengths are accounted for appropriately.
    """
    tree_1 = make_tree("(a:1,(b:1,(c:1,(d:1,e:1):1):1):1)")
    tree_2 = make_tree("(d:0.1,(f:0.1,(a:0.1,b:0.1):0.1):0.1)")

    expected_one_branch = make_tree("((f,a),(b,(c,(d,e))))")
    expected_depth = make_tree("((f,(a,b)),(c,(d,e)))")

    scs_test(
        [tree_1, tree_2],
        expected_one_branch,
        pcg_weighting="one",
        contract_edges=False,
    )
    scs_test(
        [tree_1, tree_2],
        expected_depth,
        pcg_weighting="depth",
        contract_edges=False,
    )
    scs_test(
        [tree_1, tree_2],
        expected_one_branch,
        pcg_weighting="branch",
        contract_edges=False,
    )


def test_bootstrap() -> None:
    """
    Test that bootstrap pcg weighting works appropriately
    (branch equivalent when no lengths).
    """
    tree_1 = make_tree("(a,(b,(c,(d,e)100)100)100)")
    tree_2 = make_tree("(a,(b,(d,(c,e)45)100)100)100")
    tree_3 = make_tree("(a,(b,(d,(c,e)50)100)100)100")

    expected_one_depth = make_tree("(a,(b,(d,(c,e))))")
    expected_bootstrap = make_tree("(a,(b,(c,(d,e))))")

    scs_test(
        [tree_1, tree_2, tree_3],
        expected_one_depth,
        pcg_weighting="one",
        contract_edges=False,
    )
    scs_test(
        [tree_1, tree_2, tree_3],
        expected_one_depth,
        pcg_weighting="depth",
        contract_edges=False,
    )
    scs_test(
        [tree_1, tree_2, tree_3],
        expected_bootstrap,
        pcg_weighting="bootstrap",
        contract_edges=False,
    )


def test_dcm_iq() -> None:
    expected = load_expected_tree_file("dcm_iq_expected.tre")
    source_trees = load_source_tree_file("dcm_iq_source.tre")

    scs_test(source_trees, expected, pcg_weighting="branch")


def test_supertriplets() -> None:
    expected = load_expected_tree_file("supertriplets_expected.tre")
    source_trees = load_source_tree_file("supertriplets_source.tre")

    scs_test(source_trees, expected, pcg_weighting="depth")


def test_not_completed() -> None:
    """
    Tests for when at least one tree is not completed.
    """

    tree_1 = make_tree("((a,b),(c,d))")
    tree_2 = make_tree("((a,b),(c,(d,e)))")

    expected = make_tree("((a,b),(c,(d,e)))")

    trees: list[Any] = [
        tree_1,
        tree_2,
        NotCompleted("ERROR", "local", "Example NotCompleted"),
    ]

    scs_test(trees, expected, weights=[1, 2, 3])
