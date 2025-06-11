from collections.abc import Callable, Collection, Iterable, Sequence
from typing import (
    Literal,
    NewType,
    TypeAlias,
)

import numpy as np
from cogent3 import make_tree
from cogent3.core.tree import PhyloNode, TreeBuilder, TreeNode
from sklearn.cluster import SpectralClustering

Taxa = NewType("Taxa", str)
PcgVertex: TypeAlias = tuple[Taxa, ...]
PcgVertexSet: TypeAlias = set[PcgVertex]
PcgEdgeMap: TypeAlias = dict[PcgVertex, PcgVertexSet]
EdgeTuple: TypeAlias = tuple[PcgVertex, PcgVertex]


def construct_supertree(
    trees: Sequence[TreeNode],
    weights: Sequence[float] | None = None,
    pcg_weighting: Literal["one", "branch", "depth", "bootstrap"] = "one",
    *,
    contract_edges: bool = True,
    random_state: np.random.RandomState | None = None,
) -> PhyloNode:
    """Spectral Cluster Supertree (SCS).

    Constructs a supertree from a collection of input trees. The supertree
    method is inspired by Min-Cut Supertree [1_], using
    spectral clustering instead of min-cut to improve efficiency.

    The set of input trees must overlap, the optional weights parameter
    allows the biasing of some trees over others.

    Parameters
    ----------
    trees : Sequence[TreeNode]
        The trees to find the supertree of.
    weights : Sequence[float] | None, optional
        The weights of the given trees, by default None.
    pcg_weighting : Literal["one", "branch", "depth", "bootstrap"], optional
        The weighting strategy to use, by default "one".
    contract_edges : bool, optional
        Whether to contract the edges of the proper cluster graph, by default True.
    random_state : np.random.RandomState, optional
        Random number generation to use, by default np.random.RandomState().

    Returns
    -------
    PhyloNode
        The generated supertree.

    References
    ----------
    .. [1] Semple, C., & Steel, M. (2000).
       A supertree method for rooted trees.
       Discrete Applied Mathematics, 105(1-3), 147-158.

    """
    if random_state is None:
        random_state = np.random.RandomState()

    if len(trees) == 0:
        msg = "There must be at least one tree to make a supertree."
        raise ValueError(msg)

    if pcg_weighting not in ("one", "branch", "depth", "bootstrap"):
        msg = f"Invalid weighting strategy selected: '{pcg_weighting}'"
        raise ValueError(msg)

    # Input trees are of equal weight if none is specified
    if weights is None:
        weights = [1.0 for _ in range(len(trees))]

    if len(trees) != len(weights):
        msg = (
            f"The number of trees ({len(trees)}) "
            f"and tree weights ({len(weights)}) must match."
        )
        raise ValueError(msg)

    if len(trees) == 1:  # If there is only one tree left, we can simply graft it on
        _denamify(trees[0])
        return make_tree(trees[0].get_newick())

    # The vertices of the proper cluster graph
    # are the names of the tips of all trees
    all_names = _get_all_tip_names(trees)

    # If there are less than or only two names, can instantly return a tree
    if len(all_names) <= 2:
        return _tip_names_to_tree(all_names)

    pcg_vertices: PcgVertexSet = {(name,) for name in all_names}

    (
        pcg_edges,
        pcg_weights,
        taxa_occurrences,
        taxa_co_occurrences,
    ) = _proper_cluster_graph_edges(
        pcg_vertices,
        trees,
        weights,
        pcg_weighting,
    )

    components = _get_graph_components(pcg_vertices, pcg_edges)

    if len(components) == 1:
        if contract_edges:
            # Modifies the proper cluster graph inplace
            _contract_proper_cluster_graph(
                pcg_vertices,
                pcg_edges,
                pcg_weights,
                taxa_occurrences,
                taxa_co_occurrences,
            )
        components = spectral_cluster_graph(pcg_vertices, pcg_weights, random_state)

    # The child trees corresponding to the components of the graph
    child_trees: list[PhyloNode] = []

    for pcg_component in components:
        component = _component_to_names_set(pcg_component)
        # Trivial case for if the size of the component is <=2
        # Simply add a tree expressing that
        if len(component) <= 2:
            child_trees.append(_tip_names_to_tree(component))
            continue

        # Otherwise, need to induce the trees on each component
        # and recursively call SCS

        # Note, inducing could possible remove trees.
        new_induced_trees, new_weights = _generate_induced_trees_with_weights(
            component,
            trees,
            weights,
        )

        # Find the supertree for the induced trees
        child_trees.append(
            construct_supertree(
                new_induced_trees,
                new_weights,
                pcg_weighting,
                contract_edges=contract_edges,
                random_state=random_state,
            ),
        )

        missing_tips = component.difference(_get_all_tip_names(new_induced_trees))

        # In this case, treat these tips as individual subtrees
        child_trees.extend(_tip_names_to_tree((x,)) for x in missing_tips)

    # Connect the child trees by making adjacent to a new root.
    return _connect_trees(child_trees)


def _denamify(tree: TreeNode) -> None:
    """Remove all non-tip names in the trees.

    Parameters
    ----------
    tree : TreeNode
        The trees to remove internal node names of.

    """
    for node in tree.iter_nontips(include_self=True):
        node.name = None


def _component_to_names_set(component: PcgVertexSet) -> set[Taxa]:
    """Convert the vertex representation to a set of names of taxa.

    Parameters
    ----------
    component : PcgVertexSet
        A component of the proper cluster graph.

    Returns
    -------
    set[Taxa]
        A set of names of taxa in the component.

    """
    names_set: set[Taxa] = set()
    for c in component:
        names_set.update(c)
    return names_set


def spectral_cluster_graph(
    vertices: PcgVertexSet,
    edge_weights: dict[EdgeTuple, float],
    random_state: np.random.RandomState,
) -> list[PcgVertexSet]:
    """Partition the taxa through spectral clustering.

    Given the proper cluster graph, perform Spectral Clustering
    to find the best partition of the vertices.

    Parameters
    ----------
    vertices : PcgVertexSet
        The vertices of the proper cluster graph.
    edge_weights : dict[EdgeTuple, float]
        The weights of the edges of the proper cluster graph.
    random_state : np.random.RandomState
        Random number generation for spectral clustering.

    Returns
    -------
    list[PcgVertexSet]
        The bipartition of taxa of the proper cluster graph.

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

    partition: list[PcgVertexSet] = [set(), set()]
    for vertex, idx in zip(vertex_list, idxs, strict=True):
        partition[idx].add(vertex)

    return partition


def _contract_proper_cluster_graph(
    vertices: PcgVertexSet,
    edges: PcgEdgeMap,
    edge_weights: dict[EdgeTuple, float],
    taxa_occurrences: dict[PcgVertex, int],
    taxa_co_occurrences: dict[EdgeTuple, int],
) -> None:
    """Contracts the proper cluster graph.

    This method operates in-place.

    Given the proper cluster graph, contract every edge where
    two taxa always appear together. i.e. the number of co-occurrences
    as a proper cluster is equal to the maximum number of times either
    taxa appears in any of the source trees.

    The vertices for the contracted edges is a tuple containing
    the taxa in the old vertices as elements (sorted).

    The weights for the parallel classes of edges formed through
    contraction are calculated by the maximum of the weights of the
    trees that support at least one of those edges.

    Parameters
    ----------
    vertices : PcgVertexSet
        The vertices of the proper cluster graph (modified in-place).
    edges : PcgEdgeMap
        The edges of the proper cluster graph (modified in-place).
    edge_weights : dict[EdgeTuple, float]
        The weights of the edges of the proper cluster graph (modified in-place).
    taxa_occurrences : dict[PcgVertex, int]
        The number of times each taxon appears in any of the input trees.
    taxa_co_occurrences : dict[EdgeTuple, int]
        The number of times two taxa appear as a proper cluster.

    """
    # Construct a new graph containing only the edges of maximal weight.
    # The components of this graph are the vertices following contraction
    max_vertices: PcgVertexSet = set()
    max_edges: PcgEdgeMap = {}
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

    # Find the new vertices in processed_contractions
    processed_contractions: list[PcgVertex] = []
    for contraction in contractions:
        processed: list[Taxa] = []
        for vertex in contraction:
            processed.extend(vertex)
        processed_contractions.append(tuple_sorted(processed))

    # Generate a mapping from the old to new vertices
    vertex_to_contraction: dict[PcgVertex, PcgVertex] = {}
    for contraction, new_vertex in zip(
        contractions,
        processed_contractions,
        strict=True,
    ):
        for vertex in contraction:
            vertex_to_contraction[vertex] = new_vertex

    # Contract the graph
    new_edge_weights: dict[EdgeTuple, list[float]] = {}
    for contraction, new_vertex in zip(
        contractions,
        processed_contractions,
        strict=True,
    ):
        # Remove the contraction from the graph
        vertices.difference_update(contraction)
        vertex_set = set(new_vertex)

        for vertex in contraction:
            for neighbour in edges[vertex]:
                # If the neighbour is a part of the contraction
                # Simply delete the edge weight (edge will be deleted later)
                old_edge = edge_tuple(vertex, neighbour)
                if vertex_set.issuperset(neighbour):
                    edge_weights.pop(old_edge, None)
                    continue

                # Otherwise we are connecting to something outside
                # of this contraction
                new_edge_pair = edge_tuple(
                    new_vertex,
                    vertex_to_contraction.get(neighbour, neighbour),
                )

                # There may be multiple edges to a vertex outside of the contraction
                if new_edge_pair not in new_edge_weights:
                    new_edge_weights[new_edge_pair] = []
                new_edge_weights[new_edge_pair].append(edge_weights[old_edge])

                # Delete the edge and edge weight with the neighbour
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
    for edge, weights in new_edge_weights.items():
        u, v = edge
        edges[u].add(v)
        edges[v].add(u)

        edge_weights[edge] = max(weights)


def _connect_trees(trees: Collection[PhyloNode]) -> PhyloNode:
    """Connect the input trees by making them adjacent to a new root.

    Parameters
    ----------
    trees : Collection[PhyloNode]
        The input trees to connect.

    Returns
    -------
    PhyloNode
        A tree connecting all the input trees.

    """
    if len(trees) == 1:
        (one,) = trees  # Unpack only tree
        return one
    tree_builder = TreeBuilder(constructor=PhyloNode).edge_from_edge  # type: ignore[reportArgumentType]
    return tree_builder(None, trees)


def _generate_induced_trees_with_weights(
    names: set[Taxa],
    trees: Sequence[TreeNode],
    weights: Sequence[float],
) -> tuple[list[TreeNode], list[float]]:
    """Induces the input trees on the set of names.

    A tree can be induced on a set by removing all leaves that
    are not in the set. More concisely, inducing gives a subtree
    only containing the elements in names.

    The results is a list of trees only expressing the given names
    and a list containing their corresponding weights.

    Parameters
    ----------
    names : set[Taxa]
        The taxa to induce the trees on.
    trees : Sequence[TreeNode]
        The trees to induce.
    weights : Sequence[float]
        The weights of the trees.

    Returns
    -------
    tuple[list[TreeNode], list[float]]
        The induced trees.
        The corresponding weights.

    """
    induced_trees: list[TreeNode] = []
    new_weights: list[float] = []

    for tree, weight in zip(trees, weights, strict=True):
        # If the tree would end up with less than two leaves,
        # there is no point inducing (no proper clusters)
        if len(names.intersection(tree.get_tip_names())) < 2:
            continue
        induced_trees.append(tree._get_sub_tree(names))  # noqa: SLF001 # type: ignore[reportArgumentType]
        induced_trees[-1].name = "root"
        new_weights.append(weight)

    return induced_trees, new_weights


def _get_graph_components(
    vertices: PcgVertexSet,
    edges: PcgEdgeMap,
) -> list[PcgVertexSet]:
    """Gather the components of a graph.

    Parameters
    ----------
    vertices : PcgVertexSet
        The vertices of the graph.
    edges : PcgEdgeMap
        The edges of the graph.

    Returns
    -------
    list[PcgVertexSet]
        A list of components of the proper cluster graph.

    """
    components: list[PcgVertexSet] = []

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
    pcg_vertices: PcgVertexSet,
    trees: Sequence[TreeNode],
    weights: Sequence[float],
    pcg_weighting: Literal["one", "branch", "depth", "bootstrap"],
) -> tuple[
    PcgEdgeMap,
    dict[EdgeTuple, float],
    dict[PcgVertex, int],
    dict[EdgeTuple, int],
]:
    """Construct a proper cluster graph for a collection of weighted trees.

    For a tree, two leaves belong to a proper cluster if the path connecting
    them does not pass through the root. Equivalently, they are part of a
    proper cluster if they are on the same side of the tree from the root.

    The proper cluster graph contains all the leaves of the tree as vertices.
    An edge connects two vertices if they belong to a proper cluster in any
    of the input trees. Each edge is weighted by the sum of the weights
    according to the given weighting strategy.

    Parameters
    ----------
    pcg_vertices : PcgVertexSet
        The vertices of the proper cluster graph.
    trees : Sequence[TreeNode]
        The trees to construct the proper cluster graph from.
    weights : Sequence[float]
        Associated weights of each of the trees.
    pcg_weighting : Literal["one", "branch", "depth", "bootstrap"]
        The weighting strategy to use.

    Returns
    -------
    tuple[
        PcgEdgeMap,
        dict[EdgeTuple, float],
        dict[PcgVertex, int],
        dict[EdgeTuple, int],
    ]
        The edges of the proper cluster graph.
        The weights of the edges of the proper cluster graph.
        The number of times each taxon appears in any of the input trees.
        The number of times two taxa appear as a proper cluster.

    """
    edges: PcgEdgeMap = {}
    edge_weights: dict[EdgeTuple, float] = {}
    taxa_occurrences: dict[PcgVertex, int] = {}
    taxa_co_occurrences: dict[EdgeTuple, int] = {}  # Number of times a taxa appears

    for vertex in pcg_vertices:
        edges[vertex] = set()
        taxa_occurrences[vertex] = 0

    if pcg_weighting not in ("one", "branch", "depth", "bootstrap"):
        msg = f"Invalid weighting strategy selected: '{pcg_weighting}'"
        raise ValueError(msg)

    if pcg_weighting == "one":
        length_function = lambda _length, _tree: 1  # noqa: E731
    elif pcg_weighting == "depth":
        length_function = lambda length, _tree: length + 1  # noqa: E731
    elif pcg_weighting == "branch":
        length_function = lambda length, tree: length + (  # noqa: E731
            1 if tree.length is None else tree.length
        )
    elif pcg_weighting == "bootstrap":
        length_function = lambda _length, tree: tree.params["support"]  # noqa: E731
    else:
        msg = f"Unexpected pcg weighting method '{pcg_weighting}'."
        raise ValueError(msg)

    for tree, weight in zip(trees, weights, strict=True):
        for side in tree:
            side_taxa = _dfs_pcg_weights(
                edges,
                edge_weights,
                taxa_co_occurrences,
                side,
                weight,
                0,
                length_function,
            )
            for taxa in side_taxa:
                taxa_occurrences[taxa] += 1

    return edges, edge_weights, taxa_occurrences, taxa_co_occurrences


def _dfs_pcg_weights(
    edges: PcgEdgeMap,
    edge_weights: dict[EdgeTuple, float],
    taxa_co_occurrences: dict[EdgeTuple, int],
    tree: PhyloNode,
    tree_weight: float,
    length: float,
    length_function: Callable[[float, PhyloNode], float],
) -> list[PcgVertex]:
    """Recursive helper to construct the PCG from the tree in a DFS fashion.

    As all pairs of that are a descendant of an internal but on opposite sides have
    the same wait, performing a DFS minimises computational cost of constructing
    the proper cluster graph.

    Parameters
    ----------
    edges : PcgEdgeMap
        The current edges of the proper cluster graph (modified in-place).
    edge_weights : dict[EdgeTuple, float]
        The weights of the edges of the proper cluster graph (modified in-place).
    taxa_co_occurrences : dict[EdgeTuple, int]
        The number of times two taxa appear as a proper cluster (modified in-place).
    tree : PhyloNode
        The tree/internal-node to construct the proper cluster graph from.
    tree_weight : float
        The associated weight of the tree.
    length : float
        The length from the internal node to the root of the tree.
    length_function : Callable[[float, PhyloNode], float]
        Function which applies the weighting strategy.

    Returns
    -------
    list[PcgVertex]
        All descendants of the current node.

    """
    if tree.is_tip():
        tip_name: Taxa = tree.name  # type: ignore[reportAssignmentType]
        return [(tip_name,)]

    length = length_function(length, tree)

    children_tips: list[list[PcgVertex]] = []
    for side in tree:
        child_tips = _dfs_pcg_weights(
            edges,
            edge_weights,
            taxa_co_occurrences,
            side,
            tree_weight,
            length,
            length_function,
        )
        children_tips.append(child_tips)

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
                        edge_weights.get(edge, 0) + length * tree_weight
                    )
                    taxa_co_occurrences[edge] = taxa_co_occurrences.get(edge, 0) + 1

    for i in range(1, len(children_tips)):
        children_tips[0].extend(children_tips[i])

    return children_tips[0]


def edge_tuple(v1: PcgVertex, v2: PcgVertex) -> EdgeTuple:
    """Generate an edge representing two taxa.

    Orders the taxa for consistent behaviour.

    Parameters
    ----------
    v1 : PcgVertex
        The first vertex.
    v2 : PcgVertex
        The second vertex.

    Returns
    -------
    EdgeTuple
        The unique edge representing these taxa.

    """
    if v1 < v2:
        return (v1, v2)
    return (v2, v1)


def tuple_sorted(iterable: Iterable[Taxa]) -> PcgVertex:
    """Generate a new vertex representing an iterable of taxa.

    Sorts the taxa then converts them into a tuple for predictable ordering.

    Parameters
    ----------
    iterable : Iterable[Taxa]
        An iterable of taxa.

    Returns
    -------
    PcgVertex
        A new vertex representing the group of taxa.

    """
    return tuple(sorted(iterable))


def _get_all_tip_names(trees: Iterable[TreeNode]) -> set[Taxa]:
    """Collect all taxa names from an iterable of trees.

    Parameters
    ----------
    trees : Iterable[TreeNode]
        The trees to collect the taxa of.

    Returns
    -------
    set[Taxa]
        A set containing the tip names of the trees.

    """
    names: set[Taxa] = set()
    for tree in trees:
        names.update(tree.get_tip_names())  # type: ignore[reportArgumentType]
    return names


def _tip_names_to_tree(tip_names: Iterable[Taxa]) -> PhyloNode:
    """Generate a rooted tree of the taxa.

    All tip names are made adjacent to a new root node.

    Parameters
    ----------
    tip_names : Iterable[Taxa]
        The names of the tips.

    Returns
    -------
    PhyloNode
        A star tree with a root connecting each of the tip names.

    """
    tree_builder = TreeBuilder(constructor=PhyloNode).create_edge  # type: ignore[reportArgumentType]
    tips = [tree_builder([], tip_name, {}) for tip_name in tip_names]
    return _connect_trees(tips)
