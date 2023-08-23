"""
Spectral Cluster Supertree - networkx implementation
"""

import math
import time
from typing import (
    Any,
    Collection,
    Dict,
    FrozenSet,
    Iterable,
    List,
    Optional,
    Sequence,
    Set,
    Tuple,
)

import numpy as np
import networkx as nx
from cogent3.core.tree import TreeBuilder, TreeNode
from sklearn.cluster import SpectralClustering


def spectral_cluster_supertree(
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

    # Input trees are of equal weight if none is specified
    if weights is None:
        weights = [1.0 for _ in range(len(trees))]

    assert len(trees) == len(weights), "trees and weights must be of same length"

    # The vertices of the proper cluster graph
    # are the names of the tips of all trees
    tip_names = _get_all_tip_names(trees)

    # If there are less than or only two names, can instantly return a tree
    if len(tip_names) <= 2:
        tree = _tip_names_to_tree(tip_names)
        return tree

    print("STARTING PCG")
    start = time.time()
    pcg = _proper_cluster_graph(trees, weights)
    print("TOOK", time.time() - start)

    print("STARTING COMP")
    start = time.time()
    components = list(nx.connected_components(pcg))
    print("TOOK", time.time() - start)

    if len(components) == 1:
        # TODO: If there the graph is connected, then need to perform spectral clustering
        # to find "best" components

        # Modifies the proper cluster graph inplace
        print("START CONTRACT")
        start = time.time()
        node_names = _contract_proper_cluster_graph(pcg, trees, weights)
        print("TOOK", time.time() - start)

        print("START CLUSTER")
        components = spectral_cluster_graph(pcg, node_names)
        print("TOOK", time.time() - start)

    # The child trees corresponding to the components of the graph
    child_trees = []

    # TODO: What if there are more than two components?
    # Previously I randomly resolved to make bifurcating.
    # I suppose now it makes more sense to have that as a
    # post-processing step?
    # Slightly frustrating since spectral clustering will
    # always generate two components.

    for component in components:
        # Trivial case for if the size of the component is <=2
        # Simply add a tree expressing that
        if len(component) <= 2:
            child_trees.append(_tip_names_to_tree(component))
            continue

        # Otherwise, need to induce the trees on each compoment
        # and recursively call SCS

        # Note, inducing could possible remove trees.
        new_induced_trees, new_weights = _generate_induced_trees_with_weights(
            component, trees, weights
        )

        # Find the supertree for the induced trees
        child_trees.append(spectral_cluster_supertree(new_induced_trees, new_weights))

        # It is possible that some tip names are missed (particularly
        # if inducing would only leave length 1). TODO: think more about when this case
        # if exhibited
        missing_tips = component.difference(_get_all_tip_names(new_induced_trees))

        # In this case, treat these tips as individual subtrees
        child_trees.extend(map(_tip_names_to_tree, missing_tips))

    # Connect the child trees by making adjacent to a new root.
    supertree = _connect_trees(child_trees)
    return supertree


def spectral_cluster_graph(pcg: nx.Graph, node_names: Dict) -> List[Set]:
    """
    Given the proper cluster graph, perform Spectral Clustering
    to find the best partition of the vertices.

    TODO: I asked earlier what to do as this always returns
    a bipartition, though the min-cut method could possibly
    return an arbitrary partition (every edge on any min-cut)
    was removed. Could it be possible to access the eigenvector
    to find vertices that don't neatly belong to either class?
    Some method with k-means to work out optimal number of classes
    (probably not because the normal case should be two classes).

    Args:
        vertices (Set): The set of vertices
        edge_weights (Dict[FrozenSet, float]): The weights of the edges

    Returns:
        List[Set]: A bipartition of the vertices
    """

    # TODO: assign labels also allows kmeans, and something else which looks less useful
    # Do they have an effect on the performance?
    sc = SpectralClustering(2, affinity="precomputed", assign_labels="kmeans")

    # # Order vertices
    # vertex_list = list(vertices)

    # TODO: Consider sparse matrix
    adjacency_matrix = nx.to_numpy_array(pcg)

    # print(f"SPARSITY: {1-np.count_nonzero(edges)/edges.size}")
    idxs = sc.fit_predict(adjacency_matrix)

    partition = [set(), set()]
    for node, idx in zip(pcg.nodes, idxs):
        if node in node_names:
            partition[idx].update(node_names[node])
        else:
            partition[idx].add(node)

    return partition


def _contract_proper_cluster_graph(
    pcg: nx.Graph,
    trees: Sequence[TreeNode],
    weights: Sequence[float],
) -> Dict:
    """
    This method operates in-place.

    Given the proper cluster graph, contract every edge of maximal
    weight (sum of the weights for the input trees).

    The vertices for the contracted edges is a frozenset containing
    the old vertices as elements.

    The weights for the parallel classes of edges formed through
    contraction are calculated by the sum of the weights of the trees
    that support at least one of those edges.

    Args:
        pcg (nx.Graph): The proper cluster graph
        trees (Sequence[TreeNode]): The input trees
        weights (Sequence[float]): The weights of the input trees
    """
    max_possible_weight = sum(weights)

    max_edges = filter(
        lambda e: math.isclose(
            max_possible_weight, pcg.edges[e]["weight"]
        ),  # Every edge of max weight
        pcg.edges,
    )

    # TODO double check this logic
    node_names: Dict[Any, Set] = {}
    names_to_node = {}
    for u, v in max_edges:
        # Contracts v into u
        new_u = names_to_node.get(u, u)
        new_v = names_to_node.get(v, v)

        if new_u == new_v:
            continue

        nx.contracted_nodes(pcg, new_u, new_v, self_loops=False, copy=False)
        if new_u not in node_names:
            node_names[new_u] = {new_u}
        if new_v in node_names:
            node_names[new_u].update(node_names.pop(new_v))
            for node in node_names[new_v]:
                names_to_node[node] = new_u
        node_names[new_u].add(new_v)
        names_to_node[new_v] = new_u

    edge_data = {}
    for edge in pcg.edges(node_names, data=True):
        u, v, data = edge
        if "contraction" not in data:
            continue
        # If contraction does appear in data, then we know that there were parllel edges here.
        # Need to re-weight
        u_names = node_names[u]
        v_names = node_names[v]
        edge_weight = 0
        for tree, tree_weight in zip(trees, weights):
            for child in tree:
                # TODO: efficiency here can be improved as we only need to find one element in common, not intersection
                if (
                    len(u_names.intersection(child.get_tip_names())) > 0
                    and len(v_names.intersection(child.get_tip_names())) > 0
                ):
                    # The tree supports the endpoints belonging to a proper cluster
                    edge_weight += tree_weight
        edge_data[(u, v)] = edge_weight
    nx.set_edge_attributes(pcg, edge_data, "weight")
    return node_names


def _connect_trees(trees: Collection[TreeNode]) -> TreeNode:
    """
    Connects the input trees by making them adjacent to a new root.

    Args:
        trees (Iterable[TreeNode]): The input trees to connect

    Returns:
        TreeNode: A tree connecting all the input trees
    """
    if len(trees) == 1:
        (one,) = trees  # Unpack only tree
        return one
    tree_builder = TreeBuilder(constructor=TreeNode).edge_from_edge
    return tree_builder(None, trees)


def _generate_induced_trees_with_weights(
    names: Set, trees: Sequence[TreeNode], weights: Sequence[float]
) -> Tuple[List[TreeNode], List[float]]:
    """
    Induces the input trees on the set of names.

    A tree can be induced on a set by removing all leaves that
    are not in the set. More concisely, inducing gives a subtree
    only containing the elements in names.

    The results is a list of trees only expressing the given names
    and a list containing their corresponding weights.

    Args:
        names (Set): The names to induce the trees on
        trees (List[TreeNode]): The original trees to be induced
        weights (List[float]): The corresponding weights of the trees

    Returns:
        Tuple[Sequence[TreeNode], Sequence[float]]: The induced trees and corresponding weights
    """
    induced_trees = []
    new_weights = []

    for tree, weight in zip(trees, weights):
        # If the tree would end up with less than two leaves,
        # there is no point inducing (no proper clusters)
        if len(names.intersection(tree.get_tip_names())) < 2:
            continue
        induced_trees.append(tree.get_sub_tree(names, ignore_missing=True))
        new_weights.append(weight)

    return induced_trees, new_weights


def _proper_cluster_graph(
    trees: Sequence[TreeNode], weights: Sequence[float]
) -> nx.Graph:
    """Constructs a proper cluster graph for a collection of weighted trees.

    For a tree, two leaves belong to a proper cluster if the path connecting
    them does not pass through the root. Equivalently, they are part of a
    proper cluster if they are on the same side of the tree from the root.

    The proper cluster graph contains all the leaves of the tree as vertices.
    An edge connects two vertices if they belong to a proper cluster in any
    of the input trees. Each edge is weighted by the sum of the weights of
    the trees for which the connected vertices are a proper cluster.

    Args:
        trees (Sequence[TreeNode]): The trees expressing the proper clusters
        weights (Sequence[float]): The weight of each tree
        names_to_nodes (Dict): The numbers representing each taxa name

    Returns:
        nx.Graph: The edges and weights of the edges
    """

    graphs = []

    for tree, weight in zip(trees, weights):
        # TODO: Should I error if more than two children?
        for side in tree:
            names = side.get_tip_names()

            g: nx.Graph = nx.complete_graph(names)
            nx.set_edge_attributes(g, weight, "weight")
            # edge_data = {
            #     e: pcg.edges[e]["weight"] + g.edges[e]["weight"]
            #     for e in pcg.edges & g.edges
            # }
            graphs.append(g)
            # print("HERE")
            # pcg = nx.compose(pcg, g)
            # print("THERE", g.nodes)
            # nx.set_edge_attributes(pcg, edge_data, "weight")

    pcg = nx.compose_all(graphs)

    edge_data = {}
    for g in graphs:
        g: nx.Graph
        for u, v, data in g.edges(data=True):
            edge_data[(u, v)] = edge_data.get((u, v), 0) + data["weight"]
    nx.set_edge_attributes(pcg, edge_data, "weight")
    return pcg


def _get_all_tip_names(trees: Iterable[TreeNode]) -> Set:
    """
    Fetch the tip names for some iterable of input trees.

    Args:
        trees (Iterable[TreeNode]): Input trees.

    Returns:
        Set: A set containing the tip names of the trees.
    """
    names = set()
    for tree in trees:
        names.update(tree.get_tip_names())
    return names


def _tip_names_to_tree(tip_names: Iterable) -> TreeNode:
    """
    Convert an iterable of tip names to a tree.
    The tip names are made adjacent to a new root.

    Args:
        tip_names (Iterable): the names of the tips.

    Returns:
        TreeNode: A star tree with a root connecting each of the tip names
    """
    tree_builder = TreeBuilder(
        constructor=TreeNode
    ).create_edge  # Incorrectly causes "type error" on input TreeNode due to type detection system
    tips = [tree_builder([], tip_name, {}) for tip_name in tip_names]
    # tree = tree_builder(
    #     tips, "root", {}
    # )  # Might need to change "root" to something else, it is likely only temporarily the root.
    return _connect_trees(tips)
