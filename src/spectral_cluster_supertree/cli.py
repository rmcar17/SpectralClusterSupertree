import click
from spectral_cluster_supertree import __version__, spectral_cluster_supertree
from cogent3 import make_tree, TreeNode


def load_trees(source_tree_file: str) -> list[TreeNode]:
    with open(source_tree_file, "r") as f:
        source_trees: list[TreeNode] = [make_tree(line.strip()) for line in f]  # type: ignore
    return source_trees


@click.command(no_args_is_help=True)
@click.version_option(__version__)
@click.option("-i", "--in-file", required=True, help="File containing source trees.")
@click.option("-o", "--out-file", required=True, help="Output file.")
def main(in_file: str, out_file: str):
    """
    Runs spectral cluster supertree over the given set of source trees.
    """
    source_trees = load_trees(in_file)
    supertree = spectral_cluster_supertree(source_trees)
    supertree.write(out_file, format="newick")


if __name__ == "__main__":
    main()
