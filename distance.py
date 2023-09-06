"""
The Generalized Robinson-Foulds Distance for Phylogenetic Trees

https://www.liebertpub.com/doi/10.1089/cmb.2021.0342
"""

from typing import Dict, Tuple
from cogent3.core.tree import TreeNode
from cogent3 import make_tree

from day_distance import com_clust, con_tree_cluster_table, make_psw, rename_trees


def grf_distance(tree_1: TreeNode, tree_2: TreeNode):
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
    ) = compute_cluster_related(tree_1, tree_2)

    assert len(tree_1_delta) == len(delta_H)
    delta_H_product = 0
    for key in delta_H:
        delta_H_product += (delta_H[key] + 1) * (tree_1_delta[key] + 1)

    numerator_1 = (
        (tree_2_vertices - cluster_intersection_cardinality) * tree_1_s_hat
        + tree_1_vertices * tree_2_s_hat
        - tree_1_vertices * sum_of_cluster_cardinalities
        - 2 * product_of_deltas
        + 2 * delta_H_product
    )


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
        vertices += 1
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
        else:
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

    tree_1 = tree_1.deepcopy()
    tree_2 = tree_2.deepcopy()
    inverse = rename_trees([tree_1, tree_2])

    # print(tree_1)
    # print(tree_2)
    # print(inverse)
    psws = list(map(make_psw, [tree_1, tree_2]))

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
    return cluster_intersection_cardinality, sum_of_cluster_cardinalities, delta_H


if __name__ == "__main__":
    print(compute_delta(make_tree("((x,y),(b,(c,d)))")))
    print(compute_s_hat(compute_delta(make_tree("((x,y),(b,(c,d)))"))))
    compute_cluster_related(
        make_tree("((a,b),((c,d),e),(f,(g,(h,i))),j,(k,l,m),n);"),
        make_tree("(((h,g,k,(a,b),l,f,m),i,j),(e,(c,d)),n);"),
    )
