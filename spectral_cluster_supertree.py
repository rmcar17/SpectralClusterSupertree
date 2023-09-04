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
        components = spectral_cluster_graph(pcg_vertices, pcg_weights)
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


def component_to_names_set(component: Set[Tuple]) -> Set:
    names_set = set()
    for c in component:
        names_set.update(c)
    return names_set


def spectral_cluster_graph(
    vertices: Set, edge_weights: Dict[Tuple, float]
) -> List[Set]:
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
    sc = SpectralClustering(2, affinity="precomputed", assign_labels="kmeans", n_jobs=1)

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
    # print("START SLOW?")
    # start = time.time()
    for i, v1 in enumerate(vertex_list):
        for j, v2 in enumerate(vertex_list):
            edges[i, j] = edge_weights.get(edge_tuple(v1, v2), 0)
    # print("SLOW? TOOK", time.time() - start)
    # print(f"SPARSITY: {1-np.count_nonzero(edges)/edges.size}")
    idxs = sc.fit_predict(edges)

    partition = [set(), set()]
    for vertex, idx in zip(vertex_list, idxs):
        partition[idx].add(vertex)

    return partition


def _contract_proper_cluster_graph(
    vertices: Set,
    edges: Dict,
    edge_weights: Dict[Tuple, float],
    max_weights: Dict[Tuple, float],
    trees: Sequence[TreeNode],
    weights: Sequence[float],
) -> None:
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
        vertices (Set): The set of vertices
        edges (Dict): A mapping of vertices to other vertices they connect
        edge_weights (Dict[FrozenSet, float]): The weight for each edge between two vertices
        trees (Sequence[TreeNode]): The input trees
        weights (Sequence[float]): The weights of the input trees
    """
    # max_possible_weight = sum(weights)

    # Construct a new graph containing only the edges of maximal weight.
    # The components of this graph are the vertices following contraction
    max_vertices = set()
    max_edges = {}
    for edge, weight in edge_weights.items():
        u, v = edge
        max_possible_weight = max(max_weights[u], max_weights[v])
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

    # Find the new vertices in processeded_contractions
    processed_contractions = []
    for contraction in contractions:
        processed = []
        for vertex in contraction:
            processed.extend(vertex)
        processed_contractions.append(tuple_sorted(processed))

    # Generate a mapping from the old to new vertices
    vertex_to_contraction = {}
    for contraction, new_vertex in zip(contractions, processed_contractions):
        for vertex in contraction:
            vertex_to_contraction[vertex] = new_vertex

    # Contract the graph
    new_edge_weights = {}
    for contraction, new_vertex in zip(contractions, processed_contractions):
        # Remove the contraction from the graph
        vertices.difference_update(contraction)
        vertex_set = set(new_vertex)

        for vertex in contraction:
            for neighbour in edges[vertex]:
                # If the neighbour is a part of the contraction
                # Simply delete the edge weight (edge will be deleted later)
                old_edge = edge_tuple(vertex, neighbour)
                if vertex_set.issuperset(neighbour):
                    if old_edge in edge_weights:
                        del edge_weights[old_edge]
                    continue

                # Otherwise we are connecting to something outside
                # of this contraction
                # new_edge_pair = frozenset(
                #     (
                #         new_vertex,
                #         vertex_to_contraction.get(neighbour, neighbour),
                #     )  # Be careful if the neighbour is in a different contraction
                # )
                new_edge_pair = edge_tuple(
                    new_vertex, vertex_to_contraction.get(neighbour, neighbour)
                )

                # There may be multiple edges to a vertex outside of the contraction
                # TODO: Perhaps don't store in list and have a seperate step later
                if new_edge_pair not in new_edge_weights:
                    new_edge_weights[new_edge_pair] = []
                new_edge_weights[new_edge_pair].append(edge_weights[old_edge])

                # Delete the edge and edge weight with the neighbout
                edges[neighbour].remove(vertex)
                del edge_weights[old_edge]
            # Handled all neighbours of the vertex,
            # can now delete the edges for this vertex
            del edges[vertex]

    # Can safely add the new vertices
    for vertex in processed_contractions:
        vertices.add(vertex)
        if vertex not in edges:  # Unnecessary if statement, but keeping consistent
            edges[vertex] = set()

    # Add the new edges to the graph
    for edge, weight in new_edge_weights.items():
        u, v = edge
        edges[u].add(v)
        edges[v].add(u)

        if len(new_edge_weights[edge]) == 1:
            edge_weights[edge] = new_edge_weights[edge][0]
        else:
            edge_weight = 0

            u = set(u)
            v = set(v)

            for tree, tree_weight in zip(trees, weights):
                for child in tree:
                    # TODO: efficiency here can be improved as we only need to find one element in common
                    if not u.isdisjoint(child.get_tip_names()) and not v.isdisjoint(
                        child.get_tip_names()
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
) -> Tuple[Dict[Tuple, Set[Tuple]], Dict[Tuple, float], Dict[Tuple, float]]:
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
    max_weights = {}  # Max weight possible for each vertex

    for vertex in pcg_vertices:
        edges[vertex] = set()
        max_weights[vertex] = 0

    # total = 0

    for tree, weight in zip(trees, weights):
        # TODO: Should I error if more than two children?
        for side in tree:
            names = side.get_tip_names()
            # start = time.time()
            # print("Len", len(names))
            for i in range(0, len(names)):
                ni = (names[i],)
                max_weights[ni] += weight
                for j in range(i):
                    nj = (names[j],)
                    edges[ni].add(nj)
                    edges[nj].add(ni)

                    edge = edge_tuple(ni, nj)
                    edge_weights[edge] = edge_weights.get(edge, 0) + weight
            # total += time.time() - start
    # print("CONSTRUCTION PART", total)
    # print(len(trees))
    return edges, edge_weights, max_weights


def edge_tuple(v1, v2):
    if v1 < v2:
        return (v1, v2)
    return (v2, v1)


def tuple_sorted(iterable: Iterable):
    return tuple(sorted(iterable))


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
