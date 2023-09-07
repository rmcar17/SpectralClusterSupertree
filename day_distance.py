"""
From Day's (1985) "Optimal Algorithms for Comparing Trees with Labeled Leaves"

From https://link.springer.com/content/pdf/10.1007/BF01908061.pdf
"""

from typing import List
from cogent3.core.tree import TreeNode
from cogent3 import make_tree


class PSW:
    def __init__(self) -> None:
        self.T = []
        self.M = 0  # number of interior vertices
        self.N = 0  # number of leaves

    def enter(self, v, w):
        self.T.append((v, w))
        if w == 0:  # It is a leaf
            self.N += 1
        else:
            self.M += 1

    def treset(self):
        self.index = 0

    def nvertex(self):
        if self.index >= len(self.T):
            return -1, 0
        result = self.T[self.index]
        self.index += 1
        return result

    def leftleaf(self):
        j = self.index - 1
        vj, wj = self.T[j]
        k = j - wj
        vk, wk = self.T[k]
        return vk

    def __str__(self) -> str:
        return str(self.T)


def make_psw(tree: TreeNode):
    T = PSW()
    weights = {}
    count = 1
    internal = "INTERNAL"
    for node in tree.postorder():
        if node.is_tip():
            weights[node] = 0
            T.enter(node.name, 0)
        else:
            weight = sum(map(lambda x: 1 + weights[x], node.children))
            T.enter(internal + str(count), weight)
            count += 1
            weights[node] = weight

    return T


class ClusterTable:
    def __init__(self, T: PSW) -> None:
        T.treset()
        self.X: List = [[None, None, None, None]]
        for _ in range(1, T.N):
            self.X.append([0, 0, None, None])

        leafcode = 0
        v, w = T.nvertex()
        while v != -1:
            if w == 0:
                leafcode = leafcode + 1
                R = leafcode
                self.X[int(v) - 1][2] = R
                v, w = T.nvertex()
            else:
                L = self.X[int(T.leftleaf()) - 1][2]
                v, w = T.nvertex()
                if w == 0:
                    loc = R
                else:
                    loc = L

                self.X[loc - 1][0] = L
                self.X[loc - 1][1] = R

    def encode(self, v):
        return self.X[int(v) - 1][2]

    def is_clust(self, L, R):
        return (self.X[L - 1][0] == L and self.X[L - 1][1] == R) or (
            self.X[R - 1][0] == L and self.X[R - 1][1] == R
        )

    def clear(self):
        for i in range(len(self.X)):
            if self.X[i][0] != 0 or self.X[i][1] != 0:
                self.X[i][3] = False

    def setsw(self, L, R):
        if self.X[L - 1][0] == L and self.X[L - 1][1] == R:
            self.X[L - 1][3] = True
        elif self.X[R - 1][0] == L and self.X[R - 1][1] == R:
            self.X[R - 1][3] = True

    def update(self):
        for i in range(len(self.X)):
            if not self.X[i][3]:
                self.X[i][0] = 0
                self.X[i][1] = 0

    def xreset(self):
        self.index = 0

    def nclus(self):
        self.index += 1
        while (
            self.index < len(self.X)
            and self.X[self.index][0] == 0
            and self.X[self.index][1] == 0
        ):
            self.index += 1

        if self.index >= len(self.X):
            return 0, 0

        return self.X[self.index][0], self.X[self.index][1]

    def number_of_clusters(self) -> int:
        length = 0
        for row in self.X:
            if row[0] != 0 and row[0] is not None:
                assert row[1] != 0 and row[1] is not None
                length += 1
        return length

    def __str__(self) -> str:
        return str(self.X)


def com_clust(psws: List[PSW]) -> ClusterTable:
    X = ClusterTable(psws[0])
    for i in range(1, len(psws)):
        Ti = psws[i]
        S = []
        X.clear()
        Ti.treset()
        v, w = Ti.nvertex()
        while v != -1:
            if w == 0:
                S.append((X.encode(v), X.encode(v), 1, 1))
            else:
                L, R, N, W = float("inf"), 0, 0, 1
                while w != 0:
                    Ls, Rs, Ns, Ws = S.pop()
                    L, R, N, W = min(L, Ls), max(R, Rs), N + Ns, W + Ws
                    w = w - Ws
                S.append((L, R, N, W))
                if N == R - L + 1:
                    X.setsw(L, R)
            v, w = Ti.nvertex()
        X.update()
    return X


def con_tree_psws(psws: List[PSW]) -> PSW:
    X = com_clust(psws)
    return con_tree_cluster_table(X, psws[0])


def con_tree_cluster_table(X: ClusterTable, T1: PSW) -> PSW:
    T = PSW()
    T1.treset()
    v, w = T1.nvertex()

    while v != -1:
        if w == 0:
            R = X.encode(v)
            T.enter(v, w)
            X.X[int(v) - 1][3] = len(T.T)
        else:
            L = X.encode(T1.leftleaf())
            if X.is_clust(L, R):
                T.enter(v, len(T.T) + 1 - X.X[int(T1.leftleaf()) - 1][3])
        v, w = T1.nvertex()

    return T


def normalise_trees(trees: List[TreeNode]):
    tree = trees[0]
    names = set(tree.get_tip_names())
    for other in trees[1:]:
        not_included = set(other.get_tip_names())
        not_included.difference_update(names)
        tree.extend(not_included)
        names.update(not_included)


def rename_trees(trees: List[TreeNode]):
    mapping = {}
    inverse = {}
    count = 1

    for tree in trees:
        for name in tree.get_tip_names():
            if name not in mapping:
                mapping[name] = count
                inverse[count] = name
                count += 1

    for tree in trees:
        tree.reassign_names(mapping)

    return inverse


if __name__ == "__main__":
    tree_1 = make_tree("((a,b),((c,d),e),(f,(g,(h,i))),j,(k,l,m),n);")
    tree_2 = make_tree("(((h,g,k,(a,b),l,f,m),i,j),(e,(c,d)),n);")

    # tree_1 = make_tree("(((a,b),d),(c,f));")
    # tree_2 = make_tree("((e,f),((a,b),d));")
    inverse = rename_trees([tree_1, tree_2])
    # count = 0
    # rename = {}
    # for name in tree.get_tip_names():
    #     rename[name] = count
    #     count += 1
    # print(tree)
    # tree.reassign_names(rename)

    # print(tree)

    print(tree_1)
    T1 = make_psw(tree_1)
    print(T1)
    X1 = ClusterTable(T1)
    print(X1)

    # print()

    print(tree_2)
    T2 = make_psw(tree_2)
    print(T2)
    # X2 = ClusterTable(T2)
    # print(X2)

    print()
    print(inverse)

    print(com_clust([T1, T2]))

    print()

    print(con_tree_psws([T1, T2]))

    print()

    print(inverse)
    # T.treset()
    # for i in range(T.M + T.N):
    #     v, w = T.nvertex()
    #     print(v, w)
    #     if "INTERNAL" in v:
    #         print("THING", T.leftleaf())
    #     # print(T.nvertex())

    # X = ClusterTable(T)
    # print(X)
    # X.xreset()
    # for i in range(len(X.X) + 5):
    #     print(X.nclus())
