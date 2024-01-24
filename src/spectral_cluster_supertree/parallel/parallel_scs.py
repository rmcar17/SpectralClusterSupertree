from multiprocessing.connection import Connection, wait
from typing import Collection, Optional, Sequence, Union
import multiprocessing as mp
import numpy as np
from spectral_cluster_supertree.parallel.parallel_util import (
    MergeResultOfTasks,
    Task,
    TaskDistributor,
    TaskResult,
)
import time

# from spectral_cluster_supertree.scs.scs import spectral_cluster_supertree
from cogent3 import make_tree, TreeNode

from spectral_cluster_supertree.scs.scs import spectral_cluster_supertree


def parallel_spectral_cluster_supertree(
    trees: Sequence[TreeNode],
    pcg_weighting: str = "one",
    normalise_pcg_weights: bool = False,
    depth_normalisation: bool = False,
    contract_edges: bool = True,
    weights: Optional[Sequence[float]] = None,
    num_workers: int = mp.cpu_count() - 1,
) -> TreeNode:
    """
    Spectral Cluster Supertree (SCS).

    Constructs a supertree from a collection of input trees. The supertree
    method is inspired by Min-Cut Supertree (Semple & Steel, 2000), using
    spectral clustering instead of min-cut to improve efficiency.

    The set of input trees must overlap, the optional weights parameter
    allows the biasing of some trees over others.

    Args:
        trees (Sequence[TreeNode]): Overlapping subtrees.
        weights (Optional[Sequence[float]]): Optional weights for the trees.

    Returns:
        TreeNode: The supertree containing all taxa in the input trees.
    """
    if num_workers == 0:
        return spectral_cluster_supertree(
            trees,
            pcg_weighting=pcg_weighting,
            normalise_pcg_weights=normalise_pcg_weights,
            depth_normalisation=depth_normalisation,
            contract_edges=contract_edges,
            weights=weights,
        )

    distributor = TaskDistributor(num_workers)

    distributor.initialise_workers()

    distributor.add_task(
        Task(
            "parallel_scs_split",
            (
                trees,
                pcg_weighting,
                normalise_pcg_weights,
                depth_normalisation,
                contract_edges,
                weights,
            ),
        )
    )

    return distributor.run()  # type: ignore


if __name__ == "__main__":
    # print(
    #     parallel_spectral_cluster_supertree(
    #         [make_tree("(a,(b,(c,(z,d))))"), make_tree("(b,(a,(d,(y,c))))")],
    #         num_workers=11,
    #     )
    # )

    with open("smo.6.modelTree.tre", "r") as f:
        trees = [make_tree(line.strip()) for line in f]

    start = time.time()
    supertree = parallel_spectral_cluster_supertree(
        trees,
        num_workers=0,
    )
    time_no_parallel = time.time() - start

    print(supertree)

    start = time.time()
    supertree = parallel_spectral_cluster_supertree(
        trees,
        num_workers=11,
    )
    time_parallel = time.time() - start

    print(supertree)
    print("No Parallel", time_no_parallel)
    print("Parallel", time_parallel)
