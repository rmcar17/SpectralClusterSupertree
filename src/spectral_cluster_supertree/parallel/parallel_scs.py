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
from cogent3.core.tree import TreeBuilder
from spectral_cluster_supertree.scs.scs import (
    _component_to_names_set,
    _connect_trees,
    _contract_proper_cluster_graph,
    _denamify,
    _generate_induced_trees_with_weights,
    _get_all_tip_names,
    _get_graph_components,
    _proper_cluster_graph_edges,
    _tip_names_to_tree,
    spectral_cluster_graph,
    spectral_cluster_supertree,
)

# from spectral_cluster_supertree.scs.scs import spectral_cluster_supertree
from cogent3 import make_tree, TreeNode


def func(x, y):
    return x + y + 1


def sum_sequence(sequence: Sequence[int]) -> int:
    if len(sequence) == 1:
        return sequence[0]
    midpoint = len(sequence) // 2
    return sum_sequence(sequence[:midpoint]) + sum_sequence(sequence[midpoint:])


def parallel_sum_sequence_split(
    sequence: Sequence[int],
) -> Union[int, MergeResultOfTasks]:
    print("Doing now")
    if len(sequence) <= 2**6:
        return sum(sequence)
    midpoint = len(sequence) // 2
    return MergeResultOfTasks(
        [
            Task(parallel_sum_sequence_split, (sequence[:midpoint],)),
            Task(parallel_sum_sequence_split, (sequence[midpoint:],)),
        ],
        parallel_sum_sequence_add,
    )


def parallel_sum_sequence_add(left: int, right: int) -> int:
    return left + right


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

    print("INIT DIST")
    distributor = TaskDistributor(num_workers)

    distributor.add_task(
        Task(
            parallel_scs_split,
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

    distributor.initialise_workers()

    return distributor.run()  # type: ignore


def parallel_scs_split(
    trees: Sequence[TreeNode],
    pcg_weighting: str = "one",
    normalise_pcg_weights: bool = False,
    depth_normalisation: bool = False,
    contract_edges: bool = True,
    weights: Optional[Sequence[float]] = None,
) -> Union[TreeNode, MergeResultOfTasks]:
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

    assert len(trees) >= 1, "there must be at least one tree"

    assert pcg_weighting in ["one", "branch", "depth"]

    # Input trees are of equal weight if none is specified
    if weights is None:
        weights = [1.0 for _ in range(len(trees))]

    assert len(trees) == len(weights), "trees and weights must be of same length"

    if len(trees) == 1:  # If there is only one tree left, we can simply graft it on
        _denamify(trees[0])
        return trees[0]

    # The vertices of the proper cluster graph
    # are the names of the tips of all trees
    all_names = _get_all_tip_names(trees)

    # If there are less than or only two names, can instantly return a tree
    if len(all_names) <= 2:
        tree = _tip_names_to_tree(all_names)
        return tree

    pcg_vertices = set((name,) for name in all_names)

    (
        pcg_edges,
        pcg_weights,
        taxa_ocurrences,
        taxa_co_occurrences,
    ) = _proper_cluster_graph_edges(
        pcg_vertices,
        trees,
        weights,
        pcg_weighting,
        normalise_pcg_weights,
        depth_normalisation,
    )

    components = _get_graph_components(pcg_vertices, pcg_edges)

    if len(components) == 1:
        if contract_edges:
            # Modifies the proper cluster graph inplace
            _contract_proper_cluster_graph(
                pcg_vertices,
                pcg_edges,
                pcg_weights,
                taxa_ocurrences,
                taxa_co_occurrences,
            )
        components = spectral_cluster_graph(
            pcg_vertices, pcg_weights, np.random.RandomState()
        )

    return MergeResultOfTasks(
        [
            Task(
                handle_component,
                (
                    component,
                    trees,
                    weights,
                    pcg_weighting,
                    normalise_pcg_weights,
                    depth_normalisation,
                    contract_edges,
                ),
            )
            for component in components
        ],
        _connect_trees_fix_tips,
        kwargs={"all_tip_names": all_names},
    )


def handle_component(
    component,
    trees,
    weights,
    pcg_weighting,
    normalise_pcg_weights,
    depth_normalisation,
    contract_edges,
) -> Union[TreeNode, MergeResultOfTasks]:
    component = _component_to_names_set(component)
    # Trivial case for if the size of the component is <=2
    # Simply add a tree expressing that
    if len(component) <= 2:
        return _tip_names_to_tree(component)

    # Otherwise, need to induce the trees on each compoment
    # and recursively call SCS

    # Note, inducing could possible remove trees.
    new_induced_trees, new_weights = _generate_induced_trees_with_weights(
        component, trees, weights
    )

    # Find the supertree for the induced trees
    child_tree = parallel_scs_split(
        new_induced_trees,
        pcg_weighting,
        normalise_pcg_weights,
        depth_normalisation,
        contract_edges,
        new_weights,
    )

    return child_tree


def _connect_trees_fix_tips(*args, **kwargs) -> TreeNode:
    """
    Connects the input trees by making them adjacent to a new root.

    Args:
        trees (Iterable[TreeNode]): The input trees to connect

    Returns:
        TreeNode: A tree connecting all the input trees
    """
    trees = list(args)

    all_tip_names = kwargs["all_tip_names"]

    for tree in trees:
        all_tip_names.difference_update(tree.get_tip_names())

    trees.extend(map(lambda x: _tip_names_to_tree((x,)), all_tip_names))

    if len(trees) == 1:
        (one,) = trees  # Unpack only tree
        return one
    tree_builder = TreeBuilder(constructor=TreeNode).edge_from_edge  # type: ignore
    return tree_builder(None, trees)


if __name__ == "__main__":
    print(
        parallel_spectral_cluster_supertree(
            [make_tree("(a,(b,(c,(z,d))))"), make_tree("(b,(a,(d,(y,c))))")],
            num_workers=11,
        )
    )
