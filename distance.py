"""
The Generalized Robinson-Foulds Distance for Phylogenetic Trees

https://www.liebertpub.com/doi/10.1089/cmb.2021.0342
"""

import time
from typing import Any, Dict, FrozenSet, Set, Tuple
from cogent3.core.tree import TreeNode
from cogent3 import make_tree
import networkx as nx
from day_distance import (
    ClusterTable,
    com_clust,
    con_tree_cluster_table,
    make_psw,
    rename_trees,
)


def rf_distance(tree_1: TreeNode, tree_2: TreeNode) -> float:
    tree_1_tips = set(tree_1.get_tip_names())
    tree_2_tips = set(tree_2.get_tip_names())

    assert tree_1_tips == tree_2_tips, (
        str(tree_1_tips.difference(tree_2_tips))
        + " "
        + str(tree_2_tips.difference(tree_1_tips))
    )

    tree_1 = tree_1.deepcopy()
    tree_2 = tree_2.deepcopy()
    inverse = rename_trees([tree_1, tree_2])

    psws = list(map(make_psw, [tree_1, tree_2]))

    cluster_tables = list(map(ClusterTable, psws))

    num_clusters = list(map(lambda x: x.number_of_clusters(), cluster_tables))
    intersection = 0

    cluster_table = cluster_tables[0]
    psw = psws[1]
    S = []
    psw.treset()
    v, w = psw.nvertex()
    while v != -1:
        if w == 0:
            S.append((cluster_table.encode(v), cluster_table.encode(v), 1, 1))
        else:
            L, R, N, W = float("inf"), 0, 0, 1
            while w != 0:
                Ls, Rs, Ns, Ws = S.pop()
                L, R, N, W = min(L, Ls), max(R, Rs), N + Ns, W + Ws
                w = w - Ws
            S.append((L, R, N, W))
            if N == R - L + 1 and cluster_table.is_clust(L, R):
                intersection += 1
        v, w = psw.nvertex()
    return sum(num_clusters) - 2 * intersection


def get_clusters_slow(tree: TreeNode) -> Set[FrozenSet[Any]]:
    clusters = set()
    for node in tree.postorder(include_self=True):
        clusters.add(frozenset(node.get_tip_names()))
    return clusters


def cluster_matching_distance(tree_1: TreeNode, tree_2: TreeNode) -> float:
    clusters_1 = get_clusters_slow(tree_1)
    clusters_2 = get_clusters_slow(tree_2)

    difference_1 = clusters_1.difference(clusters_2)
    difference_2 = clusters_2.difference(clusters_1)

    assert len(difference_1) == len(difference_2)

    graph = nx.Graph()

    graph.add_nodes_from(difference_1)
    graph.add_nodes_from(difference_2)

    for node_1 in difference_1:
        for node_2 in difference_2:
            graph.add_edge(
                node_1, node_2, weight=len(node_1.symmetric_difference(node_2))
            )

    distance = 0
    matching_edges = nx.min_weight_matching(graph)
    for edge in matching_edges:
        distance += graph.edges[edge]["weight"]

    return distance


def grf_distance_slow(tree_1: TreeNode, tree_2: TreeNode) -> float:
    # start = time.time()
    tree_1_clusters = get_clusters_slow(tree_1)
    tree_2_clusters = get_clusters_slow(tree_2)

    union_cardinality = len(tree_1_clusters.union(tree_2_clusters))

    numerator_1 = 0
    for x in tree_1_clusters:
        for y in tree_2_clusters.difference(tree_1_clusters):
            numerator_1 += len(x.symmetric_difference(y))

    numerator_2 = 0
    for x in tree_1_clusters.difference(tree_2_clusters):
        for y in tree_2_clusters:
            numerator_2 += len(x.symmetric_difference(y))

    denominator_1 = union_cardinality * len(tree_1_clusters)
    denominator_2 = union_cardinality * len(tree_2_clusters)

    # print("SLOW")
    # print(numerator_1, denominator_1)
    # print(numerator_2, denominator_2)
    # print("TIME", time.time() - start)
    # print(sorted(map(lambda x: sorted(tuple(x)), tree_1_clusters)))
    # print(sorted(map(lambda x: sorted(tuple(x)), tree_2_clusters)))
    return numerator_1 / denominator_1 + numerator_2 / denominator_2


def grf_distance(tree_1: TreeNode, tree_2: TreeNode) -> float:
    # start = time.time()
    tree_1_tips = set(tree_1.get_tip_names())
    tree_2_tips = set(tree_2.get_tip_names())

    m = len(tree_1_tips)
    n = len(tree_1_tips.union(tree_2_tips))
    k = len(tree_1_tips.intersection(tree_2_tips))

    tree_1_delta, tree_1_vertices = compute_delta_vertices(tree_1)
    tree_1_s_hat = compute_s_hat(tree_1_delta)

    tree_2_delta, tree_2_vertices = compute_delta_vertices(tree_2)
    tree_2_s_hat = compute_s_hat(tree_2_delta)

    product_of_deltas = compute_product_of_deltas(tree_1_delta, tree_2_delta)

    (
        cluster_intersection_cardinality,
        sum_of_cluster_cardinalities,
        delta_H,
        cluster_1_cardinality,
        cluster_2_cardinality,
        union_of_cluster_cardinality,
    ) = compute_cluster_related(tree_1, tree_2)

    assert len(tree_1_delta) == len(delta_H)
    assert len(tree_2_delta) == len(delta_H)
    delta_H_product_1 = 0
    delta_H_product_2 = 0
    for key in delta_H:
        delta_H_product_1 += (delta_H[key] + 1) * (tree_1_delta[key] + 1)
        delta_H_product_2 += (delta_H[key] + 1) * (tree_2_delta[key] + 1)

    numerator_1 = (
        (tree_2_vertices - cluster_intersection_cardinality) * tree_1_s_hat
        + tree_1_vertices * tree_2_s_hat
        - tree_1_vertices * sum_of_cluster_cardinalities
        - 2 * product_of_deltas
        + 2 * delta_H_product_1
    )

    denominator_1 = union_of_cluster_cardinality * cluster_1_cardinality

    numerator_2 = (
        (tree_1_vertices - cluster_intersection_cardinality) * tree_2_s_hat
        + tree_2_vertices * tree_1_s_hat
        - tree_2_vertices * sum_of_cluster_cardinalities
        - 2 * product_of_deltas
        + 2 * delta_H_product_2
    )

    denominator_2 = union_of_cluster_cardinality * cluster_2_cardinality

    # print("1 NUM", "DEN", numerator_1, denominator_1)
    # print("2 NUM", "DEN", numerator_2, denominator_2)

    # print("T2 vertices", tree_2_vertices)
    # print("clus intersect", cluster_intersection_cardinality)
    # print("T1 Depths", tree_1_delta)
    # print("T1 S Hat", tree_1_s_hat)
    # print("T1 vertices", tree_1_vertices)
    # print("T2 Depths", tree_2_delta)
    # print("T2 S Hat", tree_2_s_hat)
    # print("Sum of card", sum_of_cluster_cardinalities)
    # print("Prod of del", product_of_deltas)
    # print("H1", delta_H_product_1)
    # print("H2", delta_H_product_2)
    # print("UN", union_of_cluster_cardinality)
    # print("C1", cluster_1_cardinality)
    # print("C2", cluster_2_cardinality)
    # print("FAST")
    # print(numerator_1, denominator_1)
    # print(numerator_2, denominator_2)
    # print("TIME", time.time() - start)
    # print(cluster_1_cardinality, cluster_2_cardinality, union_of_cluster_cardinality)
    return numerator_1 / denominator_1 + numerator_2 / denominator_2


def compute_product_of_deltas(delta_1: Dict, delta_2: Dict):
    if len(delta_1) > len(delta_2):
        return compute_product_of_deltas(delta_2, delta_1)

    product_of_deltas = 0
    for key in delta_1:
        if key in delta_2:
            product_of_deltas += (delta_1[key] + 1) * (delta_2[key] + 1)

    return product_of_deltas


def compute_s_hat(delta: Dict) -> int:
    return sum(delta.values()) + len(delta)


def compute_delta_vertices(tree: TreeNode) -> Tuple[Dict, int]:
    delta = {}
    depth = 1
    vertices = 0

    child_index_stack = [0]
    curr = tree
    curr_children = curr.children
    curr_children_len = len(curr_children)

    while True:
        curr_index = child_index_stack[-1]
        if curr_index < curr_children_len:
            curr_child: TreeNode = curr_children[curr_index]
            if curr_child.children:
                child_index_stack.append(0)
                curr = curr_child
                curr_children = curr.children
                curr_children_len = len(curr_children)
                curr_index = 0
                depth += 1
            else:
                assert curr_child.is_tip()
                delta[curr_child.name] = depth
                child_index_stack[-1] += 1
                vertices += 1
        else:
            vertices += 1
            if curr is not tree:
                assert not curr.is_tip()
            if curr is tree:
                break
            curr = curr.parent
            curr_children = curr.children
            curr_children_len = len(curr_children)
            child_index_stack.pop()
            child_index_stack[-1] += 1
            depth -= 1
    return delta, vertices


def compute_cluster_related(tree_1: TreeNode, tree_2: TreeNode):
    assert set(tree_1.get_tip_names()) == set(tree_2.get_tip_names())
    single_clusters = len(set(tree_1.get_tip_names()))

    tree_1 = tree_1.deepcopy()
    tree_2 = tree_2.deepcopy()
    inverse = rename_trees([tree_1, tree_2])

    # print(tree_1)
    # print(tree_2)
    # print(inverse)
    psws = list(map(make_psw, [tree_1, tree_2]))

    cluster_tables = list(map(ClusterTable, psws))

    number_of_clusters = list(map(lambda x: x.number_of_clusters(), cluster_tables))

    number_of_clusters[0] += single_clusters
    number_of_clusters[1] += single_clusters

    union_of_clusters_cardinality = number_of_clusters[0]

    # print("NUMBER OF CLUS", number_of_clusters, tree_1, tree_2)
    cluster_table_1 = cluster_tables[0]
    psw_2 = psws[1]
    psw_2.treset()
    v, w = psw_2.nvertex()
    S = []
    # print(psw_2)
    while v != -1:
        # print("V", v, "W", w, S)
        if w == 0:
            S.append((cluster_table_1.encode(v), cluster_table_1.encode(v), 1, 1))
        else:
            L, R, N, W = float("inf"), 0, 0, 1
            while w != 0:
                Ls, Rs, Ns, Ws = S.pop()
                L, R, N, W = min(L, Ls), max(R, Rs), N + Ns, W + Ws
                w = w - Ws
            S.append((L, R, N, W))
            if N == R - L + 1 and cluster_table_1.is_clust(
                L, R
            ):  # Then we have found an identical cluster
                pass
                # print("FOUND CLUSTER", L, R)
            else:
                # print("FOUND DIFFERENT CLUSTER", L, R)
                union_of_clusters_cardinality += 1

        v, w = psw_2.nvertex()

    cluster_intersection = com_clust(psws)

    # TODO this doesn't count clusters of length 1 whereas the paper does.
    cluster_intersection_cardinality = 0
    sum_of_cluster_cardinalities = 0
    cluster_intersection.xreset()
    L, R = cluster_intersection.nclus()
    while L != 0:
        assert R != 0

        cluster_intersection_cardinality += 1
        sum_of_cluster_cardinalities += R - L + 1

        L, R = cluster_intersection.nclus()
    assert R == 0

    psw_tree = con_tree_cluster_table(cluster_intersection, psws[0])

    # print(psw_tree)

    stack = []
    delta_H = {}
    for i in range(len(psw_tree.T) - 1, -1, -1):
        name, weight = psw_tree.T[i]
        if weight != 0:
            if len(stack) > 0:
                stack[-1] -= weight
            stack.append(weight)
        else:
            delta_H[inverse[name]] = len(stack)
        # print(stack, name, weight)
        while len(stack) > 0 and stack[-1] == 0:
            del stack[-1]
        if i != 0:
            stack[-1] -= 1
    # print(delta_H)
    # print(psw_tree)
    return (
        cluster_intersection_cardinality + single_clusters,
        sum_of_cluster_cardinalities + single_clusters,
        delta_H,
        *number_of_clusters,
        union_of_clusters_cardinality,
    )


def expected_numerator(n, i):
    # return (n - i) * (3 * (n - i) + 1) / 2
    # Should actually be
    return (n * (3 * n - 4 * i - 1) + i * (2 * i + 4) - 4) // 2


def expected_denominator(n):
    return (2 * n - 1) * (2 * n - 2)


def caterpillar_expected(n, i):
    # Paper actually expects (n-i)(3(n-i)+1)/2(2n-1)(2n-2)
    # return (n - i) * (3 * (n - i) + 1) / (2 * (2 * n - 1) * (2 * n - 2))
    # Should actually be
    return expected_numerator(n, i) / expected_denominator(n)


if __name__ == "__main__":
    # print(compute_delta_vertices(make_tree("((x,y),(b,(c,d)))")))
    # print(compute_s_hat(compute_delta_vertices(make_tree("((x,y),(b,(c,d)))"))[0]))
    # print("STARTINg")
    # compute_cluster_related(
    #     make_tree("((a,b),((c,d),e),(f,(g,(h,i))),j,(k,l,m),n);"),
    #     make_tree("(((h,g,k,(a,b),l,f,m),i,j),(e,(c,d)),n);"),
    # )
    # print("DONE")
    # print(compute_delta_vertices(make_tree("((a,b),c)")))
    # print(grf_distance(make_tree("(a,(b,(c,d)))"), make_tree("(a,(b,(c,d)));")))

    n = 5  # Number of leaves
    i = 4  # Point at which edge contracted
    # print(grf_distance(make_tree("(a,(b,(c,(d,e))));"), make_tree("(a,(b,(c,d,e)));")))
    # print(expected_numerator(n, i))
    # print(expected_denominator(n))
    # print(caterpillar_expected(n, i))
    print(
        grf_distance(make_tree("(a,(b,(c,(d,e))));"), make_tree("(a,((b,e),(c,d)));")),
        grf_distance_slow(
            make_tree("(a,(b,(c,(d,e))));"), make_tree("(a,((b,e),(c,d)));")
        ),
    )

# T1 Clusters: {a}, {b}, {c}, {d}, {c,d}, {b,c,d}, {a,b,c,d}
# T2 Clusters: {a}, {b}, {c}, {d}, {c,d}, {a,b,c,d}
# x Clusters in T1, y clusters in T2 that aren't in T1
# First there are none
# Now x clusters in T1 that are not in T2, y clusters in T2
# x is {b,c,d}
# y is {a} {b} {c} {d} {c,d} {a,b,c,d}
# Symmetric Difference cardinality (|A| + |B| - 2 |A & B|)
# 4 2 2 2, 1, 1 -> 12. So numerator is correct. What about denominator?
# Should be cardinality of the union of the clusters, multiplied by card of C(t1) and C(t2) respectively
# T1 clusters have 7
# T2 clusters have 6
# Unions have 7
# Should be 49 and 42... what?
# Paper claims should  get 1/(2n-1) - so for n=4 should be 1/7
# (or in this case get 84 as a denominator... which is double what we have
# Ohhh, I used RF distance rather than GRF in the paper
# Paper actually expects (n-i)(3(n-i)+1)/2(2n-1)(2n-2)

# But it actually expected the numerator to be 7, not 12?
# It has a simplification where you sum

# What is the paper actually summing?
# for the RHS, it is {2,3,4} + {1}, {2,3,4} + {2}, {2,3,4} + {3},
# and skips (because sum only goes to n-1): {2,3,4} + {4}
# for the LHS, it is {2,3,4} + {1,2,3,4}, {2,3,4} + {3,4}
#  and skips {2,3,4} + {4}
# With the skips it should be:
# 4 + 2 + 2 + 1 + 1 = 10, which is still not the 7 the paper claims it should be?

# Then, the paper claims for the first part (Reminder in this example n=4, i=2)
# For the RHS, it should be n-i if j>=i else n-i+2. This checks out.
# For the LHS, it should be j-i if j>i else if j<i i-j. In this small example, checks out

# Then I think the sum is wrong. What was written
# Sum from j=1 to i-i of (i-j) + Sum from j=i+1 to n-1 of (j-i) + Sum from j=1 to i-i of (n-i+2) + Sum from j=i to n of (n-i)
# Should be
# Sum from j=1 to i-1 of (i-j) + Sum from j=i+1 to n-1 of (j-i) + Sum from j=1 to i-1 of (n-i+2) + Sum from j=i to n of (n-i)
# (last should make up for the skip, noting it wasn't in their earlier step, second sum in earlier step should go from 1 to n instead of 1 to n-1)

# It looks like they made a typo of i-i and put it in a solver, the expected numerator should instead actually be
# (n(3n-4i-1)+i(2i+4)-4)/2

#
