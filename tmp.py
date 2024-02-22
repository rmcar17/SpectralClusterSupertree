from spectral_cluster_supertree import spectral_cluster_supertree
from cogent3 import make_tree

tree_1 = make_tree("(a,(b,(c,(d,e))))")
tree_2 = make_tree("(d,(a,b))")

spectral_cluster_supertree([tree_1, tree_2], pcg_weighting="one", contract_edges=False)
