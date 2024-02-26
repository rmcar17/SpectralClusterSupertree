# Spectral Cluster Supertree

[![Coverage Status](https://coveralls.io/repos/github/rmcar17/SpectralClusterSupertree/badge.svg?branch=main)](https://coveralls.io/github/rmcar17/SpectralClusterSupertree?branch=main)

Spectral Cluster Supertree is a state-of-the-art algorithm for constructing rooted supertrees from collections of rooted source trees.

Spectral Cluster Supertree can be used on Newick formatted trees in Python in conjunction with [cogent3's](https://github.com/cogent3/cogent3) tree objects, or invoked from the command line.

Spectral Cluster Supertree can employ a number of weighting strategies that take into account the depths of nodes in the trees, as well as branch lengths. A user can specify weights of trees to add bias to some of the source trees if desired.

## Usage

### Python

```python
from spectral_cluster_supertree import load_source_trees, spectral_cluster_supertree

source_trees = load_source_trees("source_tree_file.tre")

supertree = spectral_cluster_supertree(source_trees, pcg_weighting="branch")

supertree.write("supertree_file.tre")
```

### CLI

In your environment which has spectral-cluster-supertree installed:

```bash
scs -i SOURCE_TREE_FILE -o SUPERTREE_FILE -p PCG_WEIGHTING_STRATEGY
```

## Weighting Strategies

### Proper Cluster Graph Weighting

Spectral Cluster Supertree recursively paritions the complete set of taxa to form a supertree. The core component of the algorithm involves partitioning the *proper cluster graph* through spectral clustering when the source trees are not consistent.

The *proper cluster graph* has the set of all taxa in the source trees as its vertices, and an edge connects two taxa if they appear together on the same side of the root in any of the source trees (such pairs of taxa are called **proper clusters**). Let $lca$ be the lowest common ancestor of a proper cluster. Each edge is weighted according to the specified strategy:

- **one** - The number of trees in which the pair of taxa appear as a proper cluster in.
- **depth** - The sum of the depths of the $lca$ of the proper cluster in all of the source trees.
- **branch** - The sum of the root to $lca$ branch lengths of the proper cluster in all of the source trees. If branch lengths are missing defaults to one (equivalent to depth). Do not use if some trees are missing branch lengths.

The **branch** weighting strategy is recommened where branch lengths are available. Otherwise, the **depth** weighting strategy is recommended over the **one** weighting strategy.

### Tree Weighting

In addition to the above, users may associate trees with weights to bias the results towards specific trees. Prior to the summation of the weights for an edge in the *proper cluster graph*, they are multiplied by the weight of the corresponding tree. The weight of each tree defaults to one if not specified.

An example is shown below, without the tree weights the alogrithm would randomly return either triple.

```python
>>> from spectral_cluster_supertree import spectral_cluster_supertree
>>> from cogent3 import make_tree
>>> tree_1 = make_tree("(a,(b,c))")
>>> tree_2 = make_tree("(c,(b,a))")
>>> print(spectral_cluster_supertree([tree_1, tree_2], weights=[1, 1.5]))
(c,(b,a));
```
