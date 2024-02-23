from typing import Literal
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
@click.option(
    "-p",
    "--pcg-weighting",
    help="Proper cluster graph weighting strategy.",
    default="branch",
    type=click.Choice(["one", "depth", "branch"], case_sensitive=False),
)
def main(in_file: str, out_file: str, pcg_weighting: Literal["one", "depth", "branch"]):
    """
    Runs spectral cluster supertree over the given set of source trees.
    """
    source_trees = load_trees(in_file)
    supertree = spectral_cluster_supertree(source_trees, pcg_weighting=pcg_weighting)
    supertree.write(out_file, format="newick")


if __name__ == "__main__":
    main()
