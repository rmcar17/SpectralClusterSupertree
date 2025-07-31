import os
from collections.abc import Sequence
from typing import Literal

import cogent3
import cogent3.app.typing as c3_types
import numpy as np
from cogent3.app.composable import define_app
from cogent3.util.misc import extend_docstring_from

from sc_supertree.load import load_trees as lts
from sc_supertree.scs import construct_supertree as cs


@define_app
@extend_docstring_from(lts)
def load_trees(
    source_tree_file: c3_types.IdentifierType | str | os.PathLike,  # type: ignore[reportInvalidTypeVarUse]
) -> list[cogent3.PhyloNode]:
    if not isinstance(source_tree_file, (str, os.PathLike)):
        msg = f"Invalid Path Type: '{type(source_tree_file)}'."
        raise TypeError(msg)
    return lts(source_tree_file)


@define_app
@extend_docstring_from(cs)
def sc_supertree(
    trees: list[cogent3.PhyloNode],
    weights: Sequence[float] | None = None,
    pcg_weighting: Literal["one", "branch", "depth", "bootstrap"] = "one",
    *,
    contract_edges: bool = True,
    random_state: np.random.RandomState | None = None,
) -> cogent3.PhyloNode:
    return cs(
        trees,
        weights,
        pcg_weighting,
        contract_edges=contract_edges,
        random_state=random_state,
    )


@define_app
def outgroup_root(
    tree: cogent3.PhyloNode,
    *,
    priority_outgroups: Sequence[str],
) -> cogent3.PhyloNode:
    """Outgroup root a tree.

    Outgroup roots the tree at the first outgroup found.

    Parameters
    ----------
    tree : cogent3.PhyloNode
        The tree to outgroup root.
    priority_outgroups : Sequence[str]
        A sequence of names to prioritise outgroup rooting at.
        Roots at the first name found at the sequence.

    Returns
    -------
    cogent3.PhyloNode
        A tree rooted at the first outgroup found.

    """
    tip_names = set(tree.get_tip_names())

    # Root at the first found name
    for name in priority_outgroups:
        if name in tip_names:
            return tree.rooted(name)

    msg = f"Tree does not contain any tip names in: {priority_outgroups}"
    raise ValueError(msg)
