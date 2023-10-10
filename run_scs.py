import sys
from typing import List
from cogent3 import make_tree
from cogent3.core.tree import TreeNode, PhyloNode
import os

# from spectral_cluster_supertree import (
#     spectral_cluster_supertree as spectral_cluster_supertree,
# )
from bl_spectral_cluster_supertree import (
    bl_spectral_cluster_supertree as spectral_cluster_supertree,
)


def parse_trees(file_path: str) -> List[PhyloNode]:
    trees = []
    with open(file_path, "r") as f:
        for line in f:
            line = line.strip()
            if len(line) > 0:
                trees.append(make_tree(line.strip()))
            # print(line.strip())
    return trees


def parse_dataset(dataset_name: str) -> List[PhyloNode]:
    path = "data/superfine/" + dataset_name + "/source trees/"
    file = os.listdir(path)[0]

    trees = parse_trees(path + file)

    return trees


def print_tree_stats(trees: List[PhyloNode]):
    all_names = set()
    for tree in trees:
        all_names.update(tree.get_tip_names())

    print(f"Trees: {len(trees)} Taxa: {len(all_names)}")


def output_supertree(supertree: TreeNode):
    print(supertree)
    # print(len(supertree.get_tip_names()))


def create_supertree(trees: List[PhyloNode]) -> None:
    # print_tree_stats(trees)
    supertree = spectral_cluster_supertree(trees)
    output_supertree(supertree)


def create_empirical_supertree(dataset_name: str) -> None:
    trees = parse_dataset(dataset_name)
    create_supertree(trees)


def create_simulated_supertree(file):
    trees = parse_trees(file)
    create_supertree(trees)


def create_simulated_supertrees(taxa: int, density: int):
    path = f"data/superfine/{taxa}-taxa/{density}/"
    for file in os.listdir(path):
        # print(file)
        if file.endswith(".source_trees"):
            create_simulated_supertree(path + file)


# if __name__ == "__main__":
#     option = sys.argv[1]  # e=Empirical, s=Simulated
#     assert option in ["e", "s"]

#     start_time = time.time()
#     if option == "e":
#         path_to_trees = sys.argv[2]
#         create_empirical_supertree(path_to_trees)
#     elif option == "s":
#         taxa = int(sys.argv[2])
#         density = int(sys.argv[3])
#         create_simulated_supertrees(taxa, density)
#     end_time = time.time()
# print(f"Completed in {end_time - start_time:.2f} seconds. {end_time-start_time}")

if __name__ == "__main__":
    input_trees = parse_trees(sys.argv[1])
    supertree = spectral_cluster_supertree(input_trees)
    print(supertree)
