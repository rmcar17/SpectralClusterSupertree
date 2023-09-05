from typing import Dict
from cogent3.core.tree import TreeNode
from cogent3 import make_tree


def grf_distance(tree_1: TreeNode, tree_2: TreeNode):
    tree_1_tips = set(tree_1.get_tip_names())
    tree_2_tips = set(tree_2.get_tip_names())

    m = len(tree_1_tips)
    n = len(tree_1_tips.union(tree_2_tips))
    k = len(tree_1_tips.intersection(tree_2_tips))

    tree_1_delta = compute_delta(tree_1)
    tree_1_s_hat = compute_s_hat(tree_1_delta)

    tree_2_delta = compute_delta(tree_2)
    tree_2_s_hat = compute_s_hat(tree_2_delta)


def compute_s_hat(delta: Dict) -> int:
    return sum(delta.values()) + len(delta)


def compute_delta(tree: TreeNode) -> Dict:
    delta = {}
    depth = 1

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
    return delta


if __name__ == "__main__":
    print(compute_delta(make_tree("((x,y),(b,(c,d)))")))
    print(compute_s_hat(compute_delta(make_tree("((x,y),(b,(c,d)))"))))
