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

