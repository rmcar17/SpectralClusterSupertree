import os
from collections.abc import Sequence
from functools import wraps
from typing import Literal

import cogent3
import numpy as np
from cogent3.app.composable import define_app

from sc_supertree.load import load_trees as lts
from sc_supertree.scs import construct_supertree as cs


@define_app
@wraps(lts)
def load_trees(source_tree_file: str | os.PathLike) -> list[cogent3.PhyloNode]:
    return lts(source_tree_file)


@define_app
@wraps(cs)
def sc_supertree(
    trees: Sequence[cogent3.TreeNode],
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
