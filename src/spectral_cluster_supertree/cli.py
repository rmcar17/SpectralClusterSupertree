from typing import Literal

import click

from cogent3 import TreeNode, make_tree

from spectral_cluster_supertree import __version__, spectral_cluster_supertree


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
@click.option(
    "--disable-contraction",
    help="Disable edge contraction (not recommended).",
    default=False,
    is_flag=True,
)
def scs(
    in_file: str,
    out_file: str,
    pcg_weighting: Literal["one", "depth", "branch"],
    disable_contraction: bool,
):
    """
    Runs spectral cluster supertree over the given set of source trees.
    """
    source_trees = load_trees(in_file)
    supertree = spectral_cluster_supertree(
        source_trees,
        pcg_weighting=pcg_weighting,
        contract_edges=not disable_contraction,
    )
    supertree.write(out_file, format="newick")
