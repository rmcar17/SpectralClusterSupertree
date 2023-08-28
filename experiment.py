import os
import subprocess
from cogent3 import make_tree
from cogent3.core.tree import TreeNode

data_path = "data/superfine/"
model_trees_end = ".model_tree"
source_trees_end = ".source_trees"


def simulated_experiment(taxa, density):
    path = data_path + f"{taxa}-taxa/{density}/"
    data_file_starts = []
    for file in os.listdir(path):
        if file.endswith(source_trees_end):
            data_file_starts.append(".".join(file.split(".")[:-1]))

    for data_file_start in data_file_starts:
        source_tree_file = path + data_file_start + source_trees_end
        model_tree_file = path + data_file_start + model_trees_end
        result_sup = subprocess.run(
            ["./run_superfine.sh", source_tree_file],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )

        result_scs = subprocess.run(
            ["./run_scs.sh", source_tree_file],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        # result = subprocess.run(["ls", "-l"], stdout=subprocess.PIPE)
        # print(result.stdout)
        # print(result.stderr)
        sup_tree = (
            make_tree(result_sup.stdout.decode("utf-8").strip())
            .bifurcating()
            .unrooted()
        )
        scs_tree = (
            make_tree(result_scs.stdout.decode("utf-8").strip())
            .bifurcating()
            .unrooted()
        )
        print("SUPERFINE:", sup_tree)
        print("SCS:", scs_tree)

        with open(model_tree_file, "r") as f:
            model = make_tree(f.read().strip()).bifurcating().unrooted()

        print(sup_tree.lin_rajan_moret(scs_tree))
        print(model.lin_rajan_moret(sup_tree))
        print(model.lin_rajan_moret(scs_tree))

        break


if __name__ == "__main__":
    simulated_experiment(100, 100)
