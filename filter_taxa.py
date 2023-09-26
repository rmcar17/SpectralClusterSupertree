from typing import Optional, Set
from cogent3 import make_tree
from cogent3.core.tree import TreeNode


def remove_taxa(tree: TreeNode, taxa_to_remove: Set[str]) -> Optional[TreeNode]:
    taxa_remaining = set(tree.get_tip_names()).difference(taxa_to_remove)
    if len(taxa_remaining) <= 2:
        return None
    return tree.get_sub_tree(taxa_remaining)


def filter_taxa(in_file: str, out_file: str, taxa_to_remove: Set[str]):
    with open(in_file + ".model_tree", "r") as f:
        model_tree = make_tree(f.read().strip())

    source_trees = []
    with open(in_file + ".source_trees", "r") as f:
        for line in f:
            source_trees.append(make_tree(line.strip()))

    model_tree = remove_taxa(model_tree, taxa_to_remove)

    for i in range(len(source_trees)):
        source_trees[i] = remove_taxa(source_trees[i], taxa_to_remove)

    with open(out_file + ".model_tree", "w") as f:
        f.write(str(model_tree))

    with open(out_file + ".source_trees", "w") as f:
        for source_tree in filter(lambda x: x is not None, source_trees):
            f.write(str(source_tree) + "\n")


if __name__ == "__main__":
    in_file = "data/superfine/500-taxa/20/sm_data.2"
    out_file = "test"
    # with open(in_file + ".model_tree") as f:
    #     print(make_tree(f.read().strip()).get_tip_names())
    to_remove = set()
    # for i in range(0, 40):
    #     to_remove.add("t" + str(i))

    filter_taxa(in_file, out_file, to_remove)
