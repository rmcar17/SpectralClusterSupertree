"""
Spectral Cluster Supertree

Link to Original Min-Cut Supertree: https://www.sciencedirect.com/science/article/pii/S0166218X0000202X
Link to Modified Min-Cut Supertree: https://www.researchgate.net/profile/Roderic-Page/publication/226374714_Modified_Mincut_Supertrees/links/0deec536b8f1011d10000000/Modified-Mincut-Supertrees.pdf

(I appreciate that Original and Modified are of the same length)

Potential dataset http://www.cs.utexas.edu/~phylo/datasets/supertrees.html
(Use wayback machine)
Relevant paper: https://academic.oup.com/sysbio/article/63/4/566/2848417?login=false

SMIDGen

https://sites.google.com/eng.ucsd.edu/datasets/dactalsuperfine

Superfine paper:
https://academic.oup.com/sysbio/article/61/2/214/1645715?login=true
Parallel Superfine paper:
https://dl.acm.org/doi/pdf/10.1145/2245276.2231992
Parallel Superfine paper (different authors): (one shared author rather)
https://www.sciencedirect.com/science/article/pii/S0167739X16300814 (it also talks about a bunch of other supertree methods)

Finally found superfine https://github.com/dtneves/SuperFine
"""

import math
import time
from typing import Collection, Dict, Iterable, List, Optional, Sequence, Set, Tuple

from cogent3.core.tree import TreeBuilder, TreeNode

from spectral_cluster_supertree import (
    _connect_trees,
    _contract_proper_cluster_graph,
    _generate_induced_trees_with_weights,
    _get_all_tip_names,
    _get_graph_components,
    _proper_cluster_graph_edges,
    _tip_names_to_tree,
)


def min_cut_supertree(
    trees: Sequence[TreeNode], weights: Optional[Sequence[float]] = None
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

    assert len(trees) >= 1, "there must be at least one tree"
    # print("CALLING ON", trees)

    # Input trees are of equal weight if none is specified
    if weights is None:
        weights = [1.0 for _ in range(len(trees))]

    assert len(trees) == len(weights), "trees and weights must be of same length"

    # The vertices of the proper cluster graph
    # are the names of the tips of all trees
    all_names = _get_all_tip_names(trees)

    # If there are less than or only two names, can instantly return a tree
    if len(all_names) <= 2:
        tree = _tip_names_to_tree(all_names)
        return tree

    pcg_vertices = set((name,) for name in all_names)

    # print("STARTING PCG")
    # start = time.time()
    pcg_edges, pcg_weights, max_weights = _proper_cluster_graph_edges(
        pcg_vertices, trees, weights
    )
    # print("TOOK", time.time() - start)

    # print("STARTING COMP")
    # start = time.time()
    components = _get_graph_components(pcg_vertices, pcg_edges)
    # print("TOOK", time.time() - start)

    if len(components) == 1:
        # TODO: If there the graph is connected, then need to perform spectral clustering
        # to find "best" components

        # Modifies the proper cluster graph inplace
        # print("START CONTRACT")
        # start = time.time()
        _contract_proper_cluster_graph(
            pcg_vertices, pcg_edges, pcg_weights, max_weights, trees, weights
        )
        # all_total[0] += time.time() - start
        # print("TOOK", time.time() - start)

        # print("START CLUSTER")
        # start = time.time()
        components = min_cut_partition(pcg_vertices, pcg_weights)
        # print("TOOK", time.time() - start)

    # The child trees corresponding to the components of the graph
    child_trees = []

    # TODO: What if there are more than two components?
    # Previously I randomly resolved to make bifurcating.
    # I suppose now it makes more sense to have that as a
    # post-processing step?
    # Slightly frustrating since spectral clustering will
    # always generate two components.

    # print("GOT COMPONENTS", components)

    for component in components:
        component = component_to_names_set(component)
        # Trivial case for if the size of the component is <=2
        # Simply add a tree expressing that
        if len(component) <= 2:
            child_trees.append(_tip_names_to_tree(component))
            continue

        # Otherwise, need to induce the trees on each compoment
        # and recursively call SCS

        # Note, inducing could possible remove trees.
        # print("BEFORE INDUCING", trees, "ON", component)
        new_induced_trees, new_weights = _generate_induced_trees_with_weights(
            component, trees, weights
        )
        # print("AFTER INDUCING", new_induced_trees)

        # Find the supertree for the induced trees
        child_trees.append(min_cut_supertree(new_induced_trees, new_weights))

        # It is possible that some tip names are missed (particularly
        # if inducing would only leave length 1). TODO: think more about when this case
        # if exhibited
        missing_tips = component.difference(_get_all_tip_names(new_induced_trees))

        # In this case, treat these tips as individual subtrees
        child_trees.extend(map(_tip_names_to_tree, missing_tips))

    # Connect the child trees by making adjacent to a new root.
    supertree = _connect_trees(child_trees)
    return supertree


def component_to_names_set(component: Set[Tuple]) -> Set:
    names_set = set()
    for c in component:
        names_set.update(c)
    return names_set


def min_cut_partition(vertices: Set, edge_weights: Dict):
    import networkx as nx

    graph = nx.Graph()
    graph.add_nodes_from(vertices)

    weighted_edges = []
    for edge in edge_weights:
        weighted_edges.append((*edge, edge_weights[edge]))
    graph.add_weighted_edges_from(weighted_edges)

    # print("GRAPH:", graph)
    # print("NODES:", graph.nodes)
    # print("EDGES:", graph.edges(data=True))

    cut_value, partition = nx.stoer_wagner(graph)
    # print(cut_value, partition)
    return list(map(set, partition))
