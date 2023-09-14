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


def make_trees(
    taxa: int, min_size: int, max_size: int, number_of_trees: int, times: int
):
    alignment = load_alignment(taxa)
    print(alignment)

    # TODO Can make this better (it currently always includes the tip at the root which isn't realistic)
    # Instead, try rooting the "quick tree" at every edge, choosing the one which shares the best split with
    # the model tree
    assert ROOT in alignment.names

    rename = {}
    for i, name in enumerate(set(alignment.names).difference([ROOT])):
        rename[name] = "t" + str(i)

    for i in range(number_of_trees):
        tree_names = create_source_tree_names(
            alignment.names, min_size, max_size, times
        )

        trees = []
        for tree_name in tree_names:
            subalignment: Union[Alignment, Any] = alignment.take_seqs(tree_name)
            tree: TreeNode = subalignment.quick_tree(calc="jc69")
            tree = tree.bifurcating().rooted_with_tip(ROOT)
            tree = tree.get_sub_tree(set(tree_name).difference([ROOT]))
            tree.reassign_names(rename)
            trees.append(tree)

        with open(
            f"simulated_trees/simulated_tree_data/{taxa}_taxa/{i}.source_trees", "w"
        ) as f:
            for tree in trees:
                f.write(str(tree) + "\n")

        zhu_tree = cogent3.load_tree(
            "simulated_trees/sim_align/zhu-tree-rerooted-resolved-molclock.nwk"
        )
        model_tree = zhu_tree.get_sub_tree(set(alignment.names).difference([ROOT]))
        model_tree.reassign_names(rename)

        with open(
            f"simulated_trees/simulated_tree_data/{taxa}_taxa/{i}.model_tree", "w"
        ) as f:
            f.write(str(model_tree) + "\n")


def create_source_tree_names(
    names: List, min_size: int, max_size: int, times: int
) -> List[List]:
    trees = []
    for _ in range(times):
        name_set = set([ROOT])

        names_copy = list(names)
        names_copy.remove(ROOT)
        while len(name_set) != len(names_copy) + 1:
            subtree = random.sample(names_copy, random.randint(min_size, max_size))
            name_set.update(subtree)
            subtree.append(ROOT)
            trees.append(subtree)

    return trees


if __name__ == "__main__":
    make_trees(100, 5, 20, 10, 5)
