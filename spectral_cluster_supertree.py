import math
from typing import (
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
    pcg_vertices = _get_all_tip_names(trees)

    # If there are less than or only two names, can instantly return a tree
    if len(pcg_vertices) <= 2:
        tree = _tip_names_to_tree(pcg_vertices)
        return tree

    # TODO: Construct proper cluster graph... perform spectral clustering
    pcg_edges, pcg_weights = _proper_cluster_graph_edges(pcg_vertices, trees, weights)

    components = _get_graph_components(pcg_vertices, pcg_edges)

    if len(components) == 1:
        # TODO: If there the graph is connected, then need to perform spectral clustering
        # to find "best" components

        # Modifies the proper cluster graph inplace
        _contract_proper_cluster_graph(
            pcg_vertices, pcg_edges, pcg_weights, trees, weights
        )

        components = spectral_cluster_graph(pcg_vertices, pcg_weights)

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


def spectral_cluster_graph(
    vertices: Set, edge_weights: Dict[FrozenSet, float]
) -> List[Set]:
    # TODO: assign labels also allows kmeans, and something else which looks less useful
    # Do they have an effect on the performance?
    sc = SpectralClustering(2, affinity="precomputed", assign_labels="discretize")

    # Order vertices
    vertex_list = list(vertices)

    # TODO: Previously I restricted the of the weights of the edges
    # to be ints, and I used a numpy array with dtype=np.int8 to save
    # memory. Should I somehow choose whether to use a sparse matrix
    # or not. Should I move back from now float to int? Is there some
    # kind of automatic selection that can be performed?
    edges = np.zeros((len(vertex_list), len(vertex_list)))

    # TODO: This is horridly inefficient. Generate mapping from vertices to indices
    # and iterate over edges instead as this will likely be semi-sparse
    for i, v1 in enumerate(vertex_list):
        for j, v2 in enumerate(vertex_list):
            edges[i, j] = edge_weights.get((frozenset((v1, v2))), 0)

    idxs = sc.fit_predict(edges)

    partition = [set(), set()]
    for vertex, idx in zip(vertex_list, idxs):
        if isinstance(vertex, frozenset):
            partition[idx].update(vertex)
        else:
            partition[idx].add(vertex)

    return partition


def _contract_proper_cluster_graph(
    vertices: Set,
    edges: Dict,
    edge_weights: Dict[FrozenSet, float],
    trees: Sequence[TreeNode],
    weights: Sequence[float],
) -> None:
    max_possible_weight = sum(weights)

    # Construct a new graph containing only the edges of maximal weight.
    # The components of this graph are the vertices following contraction
    max_vertices = set()
    max_edges = {}
    for edge, weight in edge_weights.items():
        if math.isclose(weight, max_possible_weight):
            # Add the connecting vertices to the graph
            for v in edge:
                max_vertices.add(v)
                if v not in max_edges:
                    max_edges[v] = set()
            u, v = edge
            max_edges[u].add(v)
            max_edges[v].add(u)

    # The components of the new graph are the new vertices after contraction
    contractions = _get_graph_components(max_vertices, max_edges)

    # Generate a mapping from the old to new vertices
    vertex_to_contraction = {}
    for contraction in contractions:
        for vertex in contraction:
            vertex_to_contraction[vertex] = frozenset(contraction)

    # Contract the graph
    new_edge_weights = {}
    for contraction in contractions:
        # Remove the contraction from the graph
        vertices.difference_update(contraction)
        new_vertex = frozenset(contraction)

        for vertex in contraction:
            for neighbour in edges[vertex]:
                # If the neighbour is a part of the contraction
                # Simply delete the edge weight (edge will be deleted later)
                if neighbour in new_vertex:
                    e = frozenset((vertex, neighbour))
                    if e in edge_weights:
                        del edge_weights[e]
                    continue

                # Otherwise we are connecting to something outside
                # of this contraction
                new_edge_pair = frozenset(
                    (
                        new_vertex,
                        vertex_to_contraction.get(neighbour, neighbour),
                    )  # Be careful if the neighbour is in a different contraction
                )

                # There may be multiple edges to a vertex outside of the contraction
                # TODO: Perhaps don't store in list and have a seperate step later
                if new_edge_pair not in new_edge_weights:
                    new_edge_weights[new_edge_pair] = []
                new_edge_weights[new_edge_pair].append(
                    edge_weights[frozenset((vertex, neighbour))]
                )

                # Delete the edge and edge weight with the neighbout
                edges[neighbour].remove(vertex)
                del edge_weights[frozenset((vertex, neighbour))]
            # Handled all neighbours of the vertex,
            # can now delete the edges for this vertex
            del edges[vertex]

    # Can safely add the new vertices
    for contraction in contractions:
        vertices.add(contraction)
        if contraction not in edges:  # Unnecessary if statement, but keeping consistent
            edges[contraction] = set()

    # Add the new edges to the graph
    for edge, weight in new_edge_weights.items():
        u, v = edge
        edges[u].add(v)
        edges[v].add(u)

        if len(new_edge_weights[edge]) == 1:
            edge_weights[edge] = new_edge_weights[edge][0]
        else:
            edge_weight = 0

            # Make sure if an edge is not contracted, it behaves as if it were a set
            # TODO: consider making all vertices a set so consistency doesn't cause issues
            # This may cause lookup annoyances however so perhaps not.
            if not isinstance(u, frozenset):
                u = frozenset((u,))
            if not isinstance(v, frozenset):
                v = frozenset((v,))

            for tree, tree_weight in zip(trees, weights):
                for child in tree:
                    # TODO: Is this precisely equivalent to the paper (re: parallel edges... particularly when both endpoints are contractions). Think it works out the same..
                    # TODO: efficiency here can be improved as we only need to find one element in common
                    if (
                        len(u.intersection(child.get_tip_names())) > 0
                        and len(v.intersection(child.get_tip_names())) > 0
                    ):
                        # The tree supports the endpoints belonging to a proper cluster
                        edge_weight += tree_weight

            edge_weights[edge] = edge_weight


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


def _get_graph_components(vertices: Set, edges: Dict) -> List[Set]:
    """
    Given a graph expressed as a set of vertices and a dictionary of
    edges (mapping vertices to sets of other vertices), find the
    components of the graph.

    Args:
        vertices (Set): The set of edges.
        edges (Dict): A mapping of vertices to sets of other vertices.

    Returns:
        List[Set]: A list of sets of vertices, each element a component.
    """
    components = []

    unexplored = vertices.copy()
    while unexplored:
        frontier = [unexplored.pop()]
        component = set(frontier)
        while frontier:
            current_vertex = frontier.pop()
            for neighbour in edges[current_vertex]:
                if neighbour not in component:
                    frontier.append(neighbour)
                    component.add(neighbour)
        components.append(component)
        unexplored.difference_update(component)

    return components


def _proper_cluster_graph_edges(
    pcg_vertices: Set, trees: Sequence[TreeNode], weights: Sequence[float]
) -> Tuple[Dict, Dict[FrozenSet, float]]:
    """Constructs a proper cluster graph for a collection of weighted trees.

    For a tree, two leaves belong to a proper cluster if the path connecting
    them does not pass through the root. Equivalently, they are part of a
    proper cluster if they are on the same side of the tree from the root.

    The proper cluster graph contains all the leaves of the tree as vertices.
    An edge connects two vertices if they belong to a proper cluster in any
    of the input trees. Each edge is weighted by the sum of the weights of
    the trees for which the connected vertices are a proper cluster.

    Args:
        pcg_vertices (Set): The names of all leaves in the input trees
        trees (Sequence[TreeNode]): The trees expressing the proper clusters
        weights (Sequence[float]): The weight of each tree

    Returns:
        Tuple[Dict, Dict[FrozenSet, float]]: The edges and weights of the edges
    """
    edges = {}
    edge_weights = {}

    for name in pcg_vertices:
        edges[name] = set()

    for tree, weight in zip(trees, weights):
        # TODO: Should I error if more than two children?
        for side in tree:
            names = side.get_tip_names()
            for i in range(1, len(names)):
                for j in range(i):
                    edges[names[i]].add(names[j])
                    edges[names[j]].add(names[i])

                    edge = frozenset((names[i], names[j]))
                    edge_weights[edge] = edge_weights.get(edge, 0) + weight

    return edges, edge_weights


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
