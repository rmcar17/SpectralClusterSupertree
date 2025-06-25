import os
from collections.abc import Sequence
from typing import Literal

import pytest
from cogent3 import TreeNode, get_app, make_tree
from cogent3.app.composable import NotCompleted
from helpers import TEST_DATA_DIR, load_expected_tree_file, load_source_tree_file


def scs_test_app(
    in_trees: Sequence[TreeNode],
    expected: TreeNode,
    weights: Sequence[float] | None = None,
    pcg_weighting: Literal["one", "branch", "depth", "bootstrap"] = "one",
    *,
    contract_edges: bool = True,
) -> None:
    app = get_app(
        "sc_supertree",
        weights=weights,
        pcg_weighting=pcg_weighting,
        contract_edges=contract_edges,
    )

    result = app(in_trees).sorted()
    expected = expected.sorted()
    assert result.same_shape(expected), f"{result} != {expected}"


def scs_test_pipeline(
    in_trees_file: str | os.PathLike,
    expected: TreeNode,
    weights: Sequence[float] | None = None,
    pcg_weighting: Literal["one", "branch", "depth", "bootstrap"] = "one",
    *,
    contract_edges: bool = True,
) -> None:
    lt = get_app("load_trees")
    scs = get_app(
        "sc_supertree",
        weights=weights,
        pcg_weighting=pcg_weighting,
        contract_edges=contract_edges,
    )
    app = lt + scs

    result = app(in_trees_file).sorted()
    expected = expected.sorted()
    assert result.same_shape(expected), f"{result} != {expected}"


def test_agreeable() -> None:
    """
    Tests for when spectral clustering is not required.
    """

    tree_1 = make_tree("((a,b),(c,d))")
    tree_2 = make_tree("((a,b),(c,(d,e)))")

    expected = make_tree("((a,b),(c,(d,e)))")

    scs_test_app([tree_1, tree_2], expected)

    tree_1 = make_tree("(((a,b),(c,d)),(z,(x,y)))")
    tree_2 = make_tree("((a,((f,g),b)),(c,(d,e)))")

    expected = make_tree("(((a,(b,(f,g))),(c,(d,e))),((x,y),z))")

    scs_test_app([tree_1, tree_2], expected)


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

    scs_test_app(source_trees, model_tree)
    scs_test_pipeline(TEST_DATA_DIR / source_tree_file, model_tree)


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
    scs_test_app([tree_1, tree_2, tree_3], expected)


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
    scs_test_app(
        [tree_1, tree_2, tree_3, tree_4, tree_5, tree_6, tree_7, tree_8],
        expected,
    )


def test_simple_contration() -> None:
    """
    Small problems where contraction is required
    """
    tree_1 = make_tree("(((a,b),c),(d,e))")
    tree_2 = make_tree("((a,b),(c,d))")

    expected = make_tree("(((a,b),c),(d,e))")
    scs_test_app([tree_1, tree_2], expected)


def test_size_two_trees() -> None:
    """
    Trees with only two leaves. No information gained but a special case nonetheless.
    """
    tree_1 = make_tree("(a,b)")
    tree_2 = make_tree("(b,c)")
    tree_3 = make_tree("(c,d)")
    tree_4 = make_tree("(b,a)")

    expected = make_tree("(a,b,c,d)")
    scs_test_app([tree_1, tree_2, tree_3], expected)

    scs_test_app([tree_1], tree_1)
    scs_test_app([tree_1, tree_1], tree_1)
    scs_test_app([tree_1, tree_4], tree_1)


def test_simple_weights() -> None:
    """
    Tests tree weighting works as expected
    """
    tree_1 = make_tree("(a,(b,c))")
    tree_2 = make_tree("(c,(a,b))")

    scs_test_app([tree_1, tree_2], tree_1, weights=[2, 1])
    scs_test_app([tree_1, tree_2], tree_1, weights=[1.001, 1])

    scs_test_app([tree_1, tree_2], tree_2, weights=[1, 2])
    scs_test_app([tree_1, tree_2], tree_2, weights=[1, 1.001])


def test_depth_weighting() -> None:
    """
    Test that depth pcg weighting works appropriately
    (branch equivalent when no lengths).
    """
    tree_1 = make_tree("(a,(b,(c,(d,e))))")
    tree_2 = make_tree("(d,(f,(a,b)))")

    expected_one = make_tree("((f,a),(b,(c,(d,e))))")
    expected_depth = make_tree("((f,(a,b)),(c,(d,e)))")

    scs_test_app(
        [tree_1, tree_2],
        expected_one,
        pcg_weighting="one",
        contract_edges=False,
    )
    scs_test_app(
        [tree_1, tree_2],
        expected_depth,
        pcg_weighting="depth",
        contract_edges=False,
    )
    scs_test_app(
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

    scs_test_app(
        [tree_1, tree_2],
        expected_one_branch,
        pcg_weighting="one",
        contract_edges=False,
    )
    scs_test_app(
        [tree_1, tree_2],
        expected_depth,
        pcg_weighting="depth",
        contract_edges=False,
    )
    scs_test_app(
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

    scs_test_app(
        [tree_1, tree_2, tree_3],
        expected_one_depth,
        pcg_weighting="one",
        contract_edges=False,
    )
    scs_test_app(
        [tree_1, tree_2, tree_3],
        expected_one_depth,
        pcg_weighting="depth",
        contract_edges=False,
    )
    scs_test_app(
        [tree_1, tree_2, tree_3],
        expected_bootstrap,
        pcg_weighting="bootstrap",
        contract_edges=False,
    )


def test_dcm_iq() -> None:
    source_trees_file = "dcm_iq_source.tre"
    expected = load_expected_tree_file("dcm_iq_expected.tre")
    source_trees = load_source_tree_file(source_trees_file)

    scs_test_app(source_trees, expected, pcg_weighting="branch")
    scs_test_pipeline(
        TEST_DATA_DIR / source_trees_file,
        expected,
        pcg_weighting="branch",
    )


def test_supertriplets() -> None:
    source_trees_file = "supertriplets_source.tre"
    expected = load_expected_tree_file("supertriplets_expected.tre")
    source_trees = load_source_tree_file(source_trees_file)

    scs_test_app(source_trees, expected, pcg_weighting="depth")
    scs_test_pipeline(
        TEST_DATA_DIR / source_trees_file,
        expected,
        pcg_weighting="depth",
    )


def _tree_equal(tree1: TreeNode, tree2: TreeNode) -> bool:
    return tree1.sorted().same_shape(tree2.sorted())


@pytest.mark.parametrize(
    "priority_outgroups",
    [("b",), ("b", "c", "d", "e", "a"), ["x", "y", "b", "c", "d"]],
)
def test_outgroup_root(priority_outgroups: Sequence[str]) -> None:
    tree = make_tree("((a,b),((c,d),(e,f)))")
    expected = make_tree("(b,(a,((c,d),(e,f))))")

    app = get_app("outgroup_root", priority_outgroups=priority_outgroups)

    got = app(tree)
    assert _tree_equal(got, expected)


@pytest.mark.parametrize(
    "priority_outgroups",
    [("c",), ["c", "b", "d", "e", "a"], ("x", "y", "c", "b", "d")],
)
def test_another_outgroup_root(priority_outgroups: Sequence[str]) -> None:
    tree = make_tree("((a,b),((c,d),(e,f)))")
    expected = make_tree("(c,(d,((e,f),(a,b))))")

    app = get_app("outgroup_root", priority_outgroups=priority_outgroups)

    got = app(tree)
    assert _tree_equal(got, expected)


@pytest.mark.parametrize(
    "priority_outgroups",
    [("x",), (), ("x", "cat")],
)
def test_missing_outgroup(priority_outgroups: Sequence[str]) -> None:
    tree = make_tree("((a,b),((c,d),(e,f)))")

    app = get_app("outgroup_root", priority_outgroups=priority_outgroups)

    got = app(tree)
    assert isinstance(got, NotCompleted)
