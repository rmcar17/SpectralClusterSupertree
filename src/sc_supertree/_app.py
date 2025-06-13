import os
from functools import wraps

import cogent3
from cogent3.app.composable import define_app

from sc_supertree.load import load_trees as lts


@define_app
@wraps(lts)
def load_trees(source_tree_file: str | os.PathLike) -> list[cogent3.PhyloNode]:
    return lts(source_tree_file)
