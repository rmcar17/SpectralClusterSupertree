from pathlib import Path

from cogent3 import PhyloNode, load_tree

from sc_supertree import load_trees

TEST_DATA_DIR = Path("tests/test_data")


def load_expected_tree_file(model_tree_file: Path | str) -> PhyloNode:
    return load_tree(TEST_DATA_DIR / model_tree_file)


def load_source_tree_file(source_tree_file: Path | str) -> list[PhyloNode]:
    return load_trees(TEST_DATA_DIR / source_tree_file)
