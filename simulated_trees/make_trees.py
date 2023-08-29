import pickle as pkl
from typing import Any, List, Union
import random
import cogent3
from cogent3.core.alignment import Alignment
from cogent3.core.tree import TreeNode
from sim_align import ROOT


def load_alignment(taxa: int) -> Alignment:
    with open("simulated_trees/aln_data/" + str(taxa) + ".pkl", "rb") as f:
        alignment = pkl.load(f)
    return alignment


def make_trees(taxa: int, min_size: int, max_size: int):
    alignment = load_alignment(taxa)
    print(alignment)

    assert ROOT in alignment.names

    rename = {}
    for i, name in enumerate(set(alignment.names).difference([ROOT])):
        rename[name] = "t" + str(i)

    tree_names = create_source_tree_names(alignment.names, min_size, max_size)

    trees = []
    for tree_name in tree_names:
        subalignment: Union[Alignment, Any] = alignment.take_seqs(tree_name)
        tree: TreeNode = subalignment.quick_tree(calc="jc69")
        tree = tree.bifurcating().rooted_with_tip(ROOT)
        tree = tree.get_sub_tree(set(tree_name).difference([ROOT]))
        tree.reassign_names(rename)
        trees.append(tree)

    with open(f"simulated_trees/simulated_tree_data/{taxa}.source_trees", "w") as f:
        for tree in trees:
            f.write(str(tree) + "\n")

    zhu_tree = cogent3.load_tree(
        "simulated_trees/sim_align/zhu-tree-rerooted-resolved-molclock.nwk"
    )
    model_tree = zhu_tree.get_sub_tree(set(alignment.names).difference([ROOT]))
    model_tree.reassign_names(rename)

    with open(f"simulated_trees/simulated_tree_data/{taxa}.model_tree", "w") as f:
        f.write(str(model_tree) + "\n")


def create_source_tree_names(names: List, min_size: int, max_size: int) -> List[List]:
    name_set = set([ROOT])

    names = list(names)
    names.remove(ROOT)

    trees = []
    while len(name_set) != len(names) + 1:
        subtree = random.sample(names, random.randint(min_size, max_size))
        name_set.update(subtree)
        subtree.append(ROOT)
        trees.append(subtree)

    return trees


if __name__ == "__main__":
    make_trees(20, 5, 10)
