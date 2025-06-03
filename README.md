# Spectral Cluster Supertree

[![PyPI Version](https://img.shields.io/pypi/v/sc-supertree)](https://pypi.org/project/sc-supertree/)
[![Python Version](https://img.shields.io/pypi/pyversions/sc-supertree)](https://pypi.org/project/sc-supertree/)
[![Code Style](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

[![CI](https://github.com/rmcar17/SpectralClusterSupertree/workflows/CI/badge.svg)](https://github.com/rmcar17/SpectralClusterSupertree/actions/workflows/ci.yml)
[![Coverage Status](https://coveralls.io/repos/github/rmcar17/SpectralClusterSupertree/badge.svg?branch=main)](https://coveralls.io/github/rmcar17/SpectralClusterSupertree?branch=main)
[![License](https://img.shields.io/github/license/rmcar17/SpectralClusterSupertree)](https://github.com/rmcar17/SpectralClusterSupertree/blob/main/LICENSE)
[![DOI](https://zenodo.org/badge/667189656.svg)](https://zenodo.org/badge/latestdoi/667189656)

Spectral Cluster Supertree is a state-of-the-art algorithm for constructing rooted supertrees from collections of rooted source trees.

Spectral Cluster Supertree can be used on Newick formatted trees in Python in conjunction with [cogent3](https://github.com/cogent3/cogent3)'s tree objects, or invoked from the command line.

Spectral Cluster Supertree can employ a number of weighting strategies that take into account the depths of nodes in the trees, as well as branch lengths. A user can specify weights of trees to add bias to some of the source trees if desired.

## Installation

```bash
pip install sc-supertree
```

## Usage

### Python

```python
from sc_supertree import load_trees, construct_supertree

source_trees = load_trees("source_tree_file.tre")

supertree = construct_supertree(source_trees, pcg_weighting="branch")

supertree.write("supertree_file.tre")
```

### CLI

In your environment which has `sc-supertree` installed:

```bash
scs -i SOURCE_TREE_FILE -o SUPERTREE_FILE -p PCG_WEIGHTING_STRATEGY
```

The ```-i``` and ```-o``` options for the input and output files are required.

The ```-p``` *proper cluster graph* weighting strategy option must be one of ```ONE|DEPTH|BRANCH|BOOTSTRAP```. It defaults to ```BRANCH``` when not provided (not recommended when some trees are missing branch lengths - see below). Tree weights are not supported through the command line.

## Weighting Strategies

### Proper Cluster Graph Weighting

Spectral Cluster Supertree recursively partitions the complete set of taxa to form a supertree. The core component of the algorithm involves partitioning the *proper cluster graph* through spectral clustering when the source trees are not consistent.

The *proper cluster graph* has the set of all taxa in the source trees as its vertices, and an edge connects two taxa if they appear together on the same side of the root in any of the source trees (such pairs of taxa are called **proper clusters**). Let $lca$ be the lowest common ancestor of a proper cluster. Each edge is weighted according to the specified strategy:

- **one** - The number of trees in which the pair of taxa appear as a proper cluster in.
- **depth** - The sum of the depths of the $lca$ of the proper cluster in all of the source trees.
- **branch** - The sum of the root to $lca$ branch lengths of the proper cluster in all of the source trees. If branch lengths are missing defaults to one (equivalent to depth). Do not use if source trees contain a mix of some trees with branch lengths and some without.
-- **bootstrap** - The sum of bootstrap values of the $lca$ nodes across trees where two taxa appear as a proper cluster.

The **branch** weighting strategy is recommended when branch lengths are available. Otherwise, the **depth** weighting strategy is recommended over the **one** weighting strategy. The **bootstrap** strategy has not yet been empirically assessed.

### Tree Weighting

In addition to the above, users may associate trees with weights to bias the results towards specific trees. Prior to the summation of the weights for an edge in the *proper cluster graph*, they are each multiplied by the weight of the corresponding tree. The weight of each tree defaults to one if not specified.

An example is shown below, without the tree weights the algorithm would randomly return either triple.

```python
>>> from sc_supertree import construct_supertree
>>> from cogent3 import make_tree
>>> tree_1 = make_tree("(a,(b,c))")
>>> tree_2 = make_tree("(c,(b,a))")
>>> print(construct_supertree([tree_1, tree_2], weights=[1, 1.5]))
(c,(b,a));
```

Tree weighting can only be used in the python implementation, not the CLI.
