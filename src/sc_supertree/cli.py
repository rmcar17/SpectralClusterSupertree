from typing import Literal

import click

from sc_supertree import __version__, construct_supertree, load_trees


@click.command(no_args_is_help=True)
@click.version_option(__version__)
@click.option("-i", "--in-file", required=True, help="File containing source trees.")
@click.option("-o", "--out-file", required=True, help="Output file.")
@click.option(
    "-p",
    "--pcg-weighting",
    help="Proper cluster graph weighting strategy.",
    default="branch",
    type=click.Choice(["one", "depth", "branch", "bootstrap"], case_sensitive=False),
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
    pcg_weighting: Literal["one", "depth", "branch", "bootstrap"],
    *,
    disable_contraction: bool,
) -> None:
    """Run spectral cluster supertree over the given set of source trees."""
    source_trees = load_trees(in_file)
    supertree = construct_supertree(
        source_trees,
        pcg_weighting=pcg_weighting,
        contract_edges=not disable_contraction,
    )
    supertree.write(out_file, format="newick")
