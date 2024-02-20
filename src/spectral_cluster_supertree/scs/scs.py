from typing import (
    Callable,
    Collection,
    Iterable,
    Optional,
    Sequence,
    NewType,
    TypeAlias,
    TypeVar,
)

import numpy as np
from cogent3.core.tree import TreeBuilder, TreeNode, PhyloNode
from sklearn.cluster import SpectralClustering


Taxa = NewType("Taxa", str)
PcgVertex: TypeAlias = tuple[Taxa, ...]
EdgeTuple: TypeAlias = tuple[PcgVertex, PcgVertex]


def spectral_cluster_supertree(
    trees: Sequence[TreeNode],
    pcg_weighting: str = "one",
    normalise_pcg_weights: bool = False,
    depth_normalisation: bool = False,
    contract_edges: bool = True,
    weights: Optional[Sequence[float]] = None,
    random_state: np.random.RandomState = np.random.RandomState(),
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

    pcg_vertices: set[PcgVertex] = set((name,) for name in all_names)

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
        components = spectral_cluster_graph(pcg_vertices, pcg_weights, random_state)

    # The child trees corresponding to the components of the graph
    child_trees: list[TreeNode] = []

    for pcg_component in components:
        component = _component_to_names_set(pcg_component)
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
        child_trees.append(
            spectral_cluster_supertree(
                new_induced_trees,
                pcg_weighting,
                normalise_pcg_weights,
                depth_normalisation,
                contract_edges,
                new_weights,
                random_state,
            )
        )

        missing_tips = component.difference(_get_all_tip_names(new_induced_trees))

        # In this case, treat these tips as individual subtrees
        child_trees.extend(map(lambda x: _tip_names_to_tree((x,)), missing_tips))

    # Connect the child trees by making adjacent to a new root.
    supertree = _connect_trees(child_trees)
    return supertree


def _denamify(tree: TreeNode):
    for node in tree.iter_nontips(include_self=True):
        node.name = None


def _component_to_names_set(component: set[PcgVertex]) -> set[Taxa]:
    names_set: set[Taxa] = set()
    for c in component:
        names_set.update(c)
    return names_set


def spectral_cluster_graph(
    vertices: set[PcgVertex],
    edge_weights: dict[EdgeTuple, float],
    random_state: np.random.RandomState,
) -> list[set[PcgVertex]]:
    """
    Given the proper cluster graph, perform Spectral Clustering
    to find the best partition of the vertices.

    Args:
        vertices (set): The set of vertices
        edge_weights (dict[Frozenset, float]): The weights of the edges

    Returns:
        list[set]: A bipartition of the vertices
    """
    sc = SpectralClustering(
        2,
        affinity="precomputed",
        assign_labels="kmeans",
        n_jobs=1,
        random_state=random_state,
    )

    # Order vertices
    vertex_list = list(vertices)

    edges = np.zeros((len(vertex_list), len(vertex_list)))

    for i, v1 in enumerate(vertex_list):
        for j, v2 in enumerate(vertex_list):
            edges[i, j] = edge_weights.get(edge_tuple(v1, v2), 0)

    idxs = sc.fit_predict(edges)

    partition: list[set[PcgVertex]] = [set(), set()]
    for vertex, idx in zip(vertex_list, idxs):
        partition[idx].add(vertex)

    return partition


def _contract_proper_cluster_graph(
    vertices: set[PcgVertex],
    edges: dict[PcgVertex, set[PcgVertex]],
    edge_weights: dict[EdgeTuple, float],
    taxa_occurrences: dict[PcgVertex, int],
    taxa_co_occurrences: dict[EdgeTuple, int],
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
        vertices (set): The set of vertices
        edges (dict): A mapping of vertices to other vertices they connect
        edge_weights (dict[Frozenset, float]): The weight for each edge between two vertices
        trees (Sequence[TreeNode]): The input trees
        weights (Sequence[float]): The weights of the input trees
    """
    # Construct a new graph containing only the edges of maximal weight.
    # The components of this graph are the vertices following contraction
    max_vertices: set[PcgVertex] = set()
    max_edges: dict[PcgVertex, set[PcgVertex]] = {}
    for pair, count in taxa_co_occurrences.items():
        u, v = pair
        max_possible_count = max(taxa_occurrences[u], taxa_occurrences[v])
        if count == max_possible_count:
            # Add the connecting vertices to the graph
            for v in pair:
                max_vertices.add(v)
                if v not in max_edges:
                    max_edges[v] = set()
            u, v = pair
            max_edges[u].add(v)
            max_edges[v].add(u)

    # The components of the new graph are the new vertices after contraction
    contractions = _get_graph_components(max_vertices, max_edges)

    # Find the new vertices in processeded_contractions
    processed_contractions: list[PcgVertex] = []
    for contraction in contractions:
        processed: list[Taxa] = []
        for vertex in contraction:
            processed.extend(vertex)
        processed_contractions.append(tuple_sorted(processed))

    # Generate a mapping from the old to new vertices
    vertex_to_contraction: dict[PcgVertex, PcgVertex] = {}
    for contraction, new_vertex in zip(contractions, processed_contractions):
        for vertex in contraction:
            vertex_to_contraction[vertex] = new_vertex

    # Contract the graph
    new_edge_weights: dict[EdgeTuple, list[float]] = {}
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
                new_edge_pair = edge_tuple(
                    new_vertex, vertex_to_contraction.get(neighbour, neighbour)
                )

                # There may be multiple edges to a vertex outside of the contraction
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
    for edge in new_edge_weights:
        u, v = edge
        edges[u].add(v)
        edges[v].add(u)

        edge_weights[edge] = max(new_edge_weights[edge])


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
    tree_builder = TreeBuilder(constructor=TreeNode).edge_from_edge  # type: ignore
    return tree_builder(None, trees)


def _generate_induced_trees_with_weights(
    names: set[Taxa], trees: Sequence[TreeNode], weights: Sequence[float]
) -> tuple[list[TreeNode], list[float]]:
    """
    Induces the input trees on the set of names.

    A tree can be induced on a set by removing all leaves that
    are not in the set. More concisely, inducing gives a subtree
    only containing the elements in names.

    The results is a list of trees only expressing the given names
    and a list containing their corresponding weights.

    Args:
        names (set): The names to induce the trees on
        trees (list[TreeNode]): The original trees to be induced
        weights (list[float]): The corresponding weights of the trees

    Returns:
        tuple[Sequence[TreeNode], Sequence[float]]: The induced trees and corresponding weights
    """
    induced_trees: list[TreeNode] = []
    new_weights: list[float] = []

    for tree, weight in zip(trees, weights):
        # If the tree would end up with less than two leaves,
        # there is no point inducing (no proper clusters)
        if len(names.intersection(tree.get_tip_names())) < 2:
            continue
        induced_trees.append(tree._get_sub_tree(names))  # type: ignore
        induced_trees[-1].name = "root"
        new_weights.append(weight)

    return induced_trees, new_weights


def _get_graph_components(
    vertices: set[PcgVertex], edges: dict[PcgVertex, set[PcgVertex]]
) -> list[set[PcgVertex]]:
    """
    Given a graph expressed as a set of vertices and a dictionary of
    edges (mapping vertices to sets of other vertices), find the
    components of the graph.

    Args:
        vertices (set): The set of edges.
        edges (dict): A mapping of vertices to sets of other vertices.

    Returns:
        list[set]: A list of sets of vertices, each element a component.
    """
    components: list[set[PcgVertex]] = []

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
    pcg_vertices: set[PcgVertex],
    trees: Sequence[TreeNode],
    weights: Sequence[float],
    pcg_weighting: str,
    normalise_pcg_weights: bool,
    depth_normalisation: bool,
) -> tuple[
    dict[PcgVertex, set[PcgVertex]],
    dict[EdgeTuple, float],
    dict[PcgVertex, int],
    dict[EdgeTuple, int],
]:
    """Constructs a proper cluster graph for a collection of weighted trees.

    For a tree, two leaves belong to a proper cluster if the path connecting
    them does not pass through the root. Equivalently, they are part of a
    proper cluster if they are on the same side of the tree from the root.

    The proper cluster graph contains all the leaves of the tree as vertices.
    An edge connects two vertices if they belong to a proper cluster in any
    of the input trees. Each edge is weighted by the sum of the weights of
    the trees for which the connected vertices are a proper cluster.

    Args:
        pcg_vertices (set): The names of all leaves in the input trees
        trees (Sequence[TreeNode]): The trees expressing the proper clusters
        weights (Sequence[float]): The weight of each tree

    Returns:
        tuple[dict, dict[Frozenset, float]]: The edges and weights of the edges
    """
    edges: dict[PcgVertex, set[PcgVertex]] = {}
    edge_weights: dict[EdgeTuple, float] = {}
    taxa_occurrences: dict[PcgVertex, int] = {}
    taxa_co_occurrences: dict[EdgeTuple, int] = {}  # Number of times a taxa appears

    for vertex in pcg_vertices:
        edges[vertex] = set()
        taxa_occurrences[vertex] = 0

    assert pcg_weighting in ["one", "branch", "depth"]

    if pcg_weighting == "one":
        length_function = lambda length, tree: 1
    elif pcg_weighting == "depth":
        length_function = lambda length, tree: length + 1
    else:
        length_function = lambda length, tree: length + (
            1 if tree.length is None else tree.length
        )

    normalise_length = 0.0
    for tree, weight in zip(trees, weights):
        depth_normalisation_factor = 1
        if depth_normalisation:
            max_length = 0
            for tip in tree.tips():
                length = 0
                while tip.parent is not None:
                    if hasattr(tip, "length"):
                        length += getattr(tip, "length")
                    else:
                        length += 1
                    tip = tip.parent
                max_length = max(max_length, length)
            depth_normalisation_factor = max_length
        for side in tree:
            side_taxa, max_sublength = dfs_pcg_weights(
                edges,
                edge_weights,
                taxa_co_occurrences,
                side,
                weight,
                0,
                length_function,
                depth_normalisation_factor,
            )
            for taxa in side_taxa:
                taxa_occurrences[taxa] += 1
            normalise_length = max(normalise_length, max_sublength)
    if normalise_pcg_weights:
        for edge in edge_weights:
            edge_weights[edge] /= normalise_length

    return edges, edge_weights, taxa_occurrences, taxa_co_occurrences


def dfs_pcg_weights(
    edges: dict[PcgVertex, set[PcgVertex]],
    edge_weights: dict[EdgeTuple, float],
    taxa_co_occurrences: dict[EdgeTuple, int],
    tree: PhyloNode,
    tree_weight: float,
    length: float,
    length_function: Callable[[float, PhyloNode], float],
    depth_normalisation_factor: int,
) -> tuple[list[PcgVertex], float]:
    if tree.is_tip():
        tip_name: Taxa = tree.name  # type: ignore
        return [(tip_name,)], 0.0

    length = length_function(length, tree)

    max_length = length
    children_tips: list[list[PcgVertex]] = []
    for side in tree:
        child_tips, normalise_length = dfs_pcg_weights(
            edges,
            edge_weights,
            taxa_co_occurrences,
            side,
            tree_weight,
            length,
            length_function,
            depth_normalisation_factor,
        )
        children_tips.append(child_tips)
        max_length = max(max_length, normalise_length)

    # Add edges to the graph
    for i in range(1, len(children_tips)):
        child_tips_1 = children_tips[i]
        for j in range(i):
            child_tips_2 = children_tips[j]

            for taxa_1 in child_tips_1:
                for taxa_2 in child_tips_2:
                    edges[taxa_1].add(taxa_2)
                    edges[taxa_2].add(taxa_1)

                    edge = edge_tuple(taxa_1, taxa_2)
                    edge_weights[edge] = (
                        edge_weights.get(edge, 0)
                        + length * tree_weight / depth_normalisation_factor
                    )
                    taxa_co_occurrences[edge] = taxa_co_occurrences.get(edge, 0) + 1

    for i in range(1, len(children_tips)):
        children_tips[0].extend(children_tips[i])

    return children_tips[0], max_length


def edge_tuple(v1: PcgVertex, v2: PcgVertex) -> EdgeTuple:
    if v1 < v2:
        return (v1, v2)
    return (v2, v1)


def tuple_sorted(iterable: Iterable[Taxa]) -> PcgVertex:
    return tuple(sorted(iterable))


def _get_all_tip_names(trees: Iterable[TreeNode]) -> set[Taxa]:
    """
    Fetch the tip names for some iterable of input trees.

    Args:
        trees (Iterable[TreeNode]): Input trees.

    Returns:
        set: A set containing the tip names of the trees.
    """
    names: set[Taxa] = set()
    for tree in trees:
        names.update(tree.get_tip_names())  # type: ignore
    return names


def _tip_names_to_tree(tip_names: Iterable[Taxa]) -> TreeNode:
    """
    Convert an iterable of tip names to a tree.
    The tip names are made adjacent to a new root.

    Args:
        tip_names (Iterable): the names of the tips.

    Returns:
        TreeNode: A star tree with a root connecting each of the tip names
    """
    tree_builder = TreeBuilder(constructor=TreeNode).create_edge  # type: ignore
    tips = [tree_builder([], tip_name, {}) for tip_name in tip_names]
    return _connect_trees(tips)
