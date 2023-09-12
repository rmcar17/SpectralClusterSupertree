from cogent3.core.tree import TreeNode
from cogent3 import make_tree


def rooted_nni(tree: TreeNode, nni_id: int, copy: bool = True) -> TreeNode:
    """
    Given a rooted phylogenetic tree and a nearest neighbour interchange identifier,
    performs the relevant nearest neighbour interchange.

    If copy is True, returns a new phylogenetic tree. Otherwise operates in-place.

    A rooted binary phylogenetic tree with n leaves has n-2 internal nodes excluding the root.
    A rooted NNI operation can be performed by, at a given internal node, swapping one of the
    children with the internal node's sibling. There are two children, so each non-root internal
    node thus corresponds to two NNI operations.

    nni_id // 2 gives the internal node in a post-order traversal of the non-root internal nodes
    to perform the NNI operation at

    If nni_id % 2 is 0, swaps the left child with the sibling. Otherwise swaps the right child.


    Args:
        tree (TreeNode): The tree to perform NNI on.
        nni_id (int): An identifier representing where to perform the NNI operation.
        copy (bool, optional): If True, operates on a copy. Otherwise in-place. Defaults to True.

    Returns:
        TreeNode: The resulting phylogenetic tree
    """
    if copy:
        tree = tree.deepcopy()

    internal_node_position = nni_id // 2
    child_index = nni_id % 2

    success = False
    position = 0
    for node in tree.postorder():
        if node.is_tip() or node.is_root():
            continue
        if position == internal_node_position:
            assert len(node.children) == 2 and node.parent is not None
            parent: TreeNode = node.parent
            assert len(parent.children) == 2
            sibling_index = 0 if parent.children[0] is not node else 1
            assert parent.children[1 - sibling_index] is node, (
                str(node) + " " + str(parent.children[1 - sibling_index])
            )
            old_sibling = parent.children[sibling_index]
            parent.children[sibling_index] = node.children[child_index]
            node.children[child_index] = old_sibling

            parent.children[sibling_index].parent = parent
            node.children[child_index].parent = node

            success = True
            break
        position += 1

    assert success, f"nni_id={nni_id} is invalid"

    return tree


def num_rooted_nni_operations(tree_size: int) -> int:
    """
    Given the size (number of leaves) of a rooted binary tree,
    returns the number of NNI operations.

    Args:
        tree_size (int): Size of a rooted binary phylogenetic tree.

    Returns:
        int: The number of NNI operations
    """
    return 2 * (tree_size - 2)


if __name__ == "__main__":
    # tree = make_tree("((a,b),(c,(d,e)));")
    # print(tree, "->", rooted_nni(tree, 5))
    tree = make_tree("((t194,t002),((t095,t117),((t190,t175),(t064,t184))));")
    for i in range(num_rooted_nni_operations(len(tree.get_tip_names()))):
        print("DOING", i)
        print(rooted_nni(tree, i, False))
