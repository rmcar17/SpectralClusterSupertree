"""
Spectral Cluster Supertree

Link to Original Min-Cut Supertree: https://pdf.sciencedirectassets.com/271602/1-s2.0-S0166218X00X00820/1-s2.0-S0166218X0000202X/main.pdf?X-Amz-Security-Token=IQoJb3JpZ2luX2VjEB0aCXVzLWVhc3QtMSJHMEUCIDx6%2Bq%2BMp49c%2B0DpWKQGgT41CHeS7uvmhrGcnT%2Fz4sUiAiEApxQbkKmRQpmX%2FAryjbpCO%2FjAeMgT3PerU%2B%2Bait5LsRgqvAUItv%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FARAFGgwwNTkwMDM1NDY4NjUiDAKOivgr26PPbbzj3SqQBWc4KFtzM4hABM7UT1TphDNvU9t6FJ4Af5L3f4sB6BfJKBaaCOjYAIzqbe1QtXJ3fQH8JEi%2F1WTsE4PKYwhEho3%2BNxaOAKtUKs%2F5QRwq67jaA9DrXAGsLnpylbrh%2FJh6kSxLjwoOEal7Id09q%2Fbu9YjFDex0J8NjnEdtcrX14QGijoErz9w9imlH8v3SFZdlomVY3HG3XY6rRUavnditiiHIS1zoKMQ2%2FkLa%2FlvUQlqdEgypifq4ca5mG%2F5BuIXDAYt04y%2B3khd6IIsvmi6IURu5zUAK57jOfx9gSCGjo%2B7EUX1wmkHoNsQEZhPVyReNwDFtSyHV0vFVZnrpHdAnY1ccUe%2FZKq3sJF59Fz4btvekrkHquBghvoroKBcdcAPT2ATO2qCKk1c0RHUSSr%2B6xBJgNNThzlwDQM0j56GZR64h%2FRyUk6MmNrGG7Hupt4r8N0pIhB1dmTBHJPYKU%2FxfrMhGvOzfwQphBt%2BDqFiULm74k3WxQrw5zc%2FofHbUWouyFUMo%2BgX7UCcXbfO9qtqJ4OFRQ%2BP7URAtthmEEPl0e8mJz7gcaeeEqj6%2BFsZ7m4LPYrVAww3x5qVQqtbxBMQPtxpwtw3rKaSdxjt6aNX39%2F9Pl9x7Pk8NS1VgHA2jLPdVghDfTUhPG8QsWFWbb6vBdXkIkeFed9O51F2qol95bNzm5sS845JExWl0SIPSm0T4fwHCCUUn%2BwCC5lJSVYrvG5o2DwW93UvPuSR3zFkxR3pLQd%2FZp2%2B%2BNLumVOR%2F3XLBQ7AKvNlko9w7W1LWChYy0QEOWPurDNAF0mD%2FAPDstzx%2Bqygik7ra2DmK%2FhArHqvr%2F1iShZzRSFIW1MJNwDEjc4alsK9JG1yWAPfH4SWZ2dSkMJX5nKYGOrEBIbVfvnec7DfEZb2h9sF15i6N%2B637BNtOuI6O0kUnYrNq38qgxPKriSHaL3Y5Oa5WqczDPtD17tO7q%2F7bRFQSuj%2FfpMOu3UvoLK2rmTw3B6RE0ZARjfaOibvJ8ZJbF8KW%2BrsgsuZPo2QMERH2gxxM%2BEnKfPWFIsnhZ%2FSJHg9vPKvihxKqaYaE41nCTMqjoX7XH0fIwjoL5FxvuGexX4dsTMHvx2Ln4nk0ECYOQjVeMgqR&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20230731T053622Z&X-Amz-SignedHeaders=host&X-Amz-Expires=300&X-Amz-Credential=ASIAQ3PHCVTYQM5ULNKL%2F20230731%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=be63c1722fca0ec03d9d663d338b01c2f26c01b56fbda122b137951653187398&hash=c311b130010be1c6d77894cdf00d9f11153a4fd4bc34f6772298868bfbe13043&host=68042c943591013ac2b2430a89b270f6af2c76d8dfd086a07176afe7c76c2c61&pii=S0166218X0000202X&tid=spdf-f88ba228-921f-40f5-81bf-9d446bbb8f7a&sid=b2a930266e4f88448779cef50668c773cc80gxrqa&type=client&tsoh=d3d3LnNjaWVuY2VkaXJlY3QuY29t&ua=07165805055450000402&rr=7ef37bba3bb955bf&cc=au
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
        weights = [1000.0 for _ in range(len(trees))]

    assert len(trees) == len(weights), "trees and weights must be of same length"

    # The vertices of the proper cluster graph
    # are the names of the tips of all trees
    pcg_vertices = _get_all_tip_names(trees)

    # If there are less than or only two names, can instantly return a tree
    if len(pcg_vertices) <= 2:
        tree = _tip_names_to_tree(pcg_vertices)
        return tree

    print("STARTING PCG")
    start = time.time()
    pcg_edges, pcg_weights = _proper_cluster_graph_edges(pcg_vertices, trees, weights)
    print("TOOK", time.time() - start)

    print("STARTING COMP")
    start = time.time()
    components = _get_graph_components(pcg_vertices, pcg_edges)
    print("TOOK", time.time() - start)

    if len(components) == 1:
        # TODO: If there the graph is connected, then need to perform spectral clustering
        # to find "best" components

        # Modifies the proper cluster graph inplace
        print("START CONTRACT")
        start = time.time()
        _contract_proper_cluster_graph(
            pcg_vertices, pcg_edges, pcg_weights, trees, weights
        )
        print("TOOK", time.time() - start)

        print("START CLUSTER")
        components = spectral_cluster_graph(pcg_vertices, pcg_weights)
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


def spectral_cluster_graph(
    vertices: Set, edge_weights: Dict[FrozenSet, float]
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
    sc = SpectralClustering(2, affinity="precomputed", assign_labels="kmeans")

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
    start = time.time()
    for i, v1 in enumerate(vertex_list):
        for j, v2 in enumerate(vertex_list):
            edges[i, j] = edge_weights.get((frozenset((v1, v2))), 0)
    # print("SLOW? TOOK", time.time() - start)
    # print(f"SPARSITY: {1-np.count_nonzero(edges)/edges.size}")
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
        c = frozenset(contraction)
        vertices.add(c)
        if c not in edges:  # Unnecessary if statement, but keeping consistent
            edges[c] = set()

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

    # total = 0

    for tree, weight in zip(trees, weights):
        # TODO: Should I error if more than two children?
        for side in tree:
            names = side.get_tip_names()
            # start = time.time()
            # print("Len", len(names))
            for i in range(1, len(names)):
                for j in range(i):
                    edges[names[i]].add(names[j])
                    edges[names[j]].add(names[i])

                    edge = frozenset((names[i], names[j]))
                    edge_weights[edge] = edge_weights.get(edge, 0) + weight
            # total += time.time() - start
    # print("CONSTRUCTION PART", total)
    # print(len(trees))
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
