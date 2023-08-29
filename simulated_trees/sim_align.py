import pickle
import sys

import cogent3
from numpy.random import default_rng

rng = default_rng(seed=1)  # setting seed will allow for reproducibility

ROOT = "G001280565"


def get_sim_params():
    """fits a ssGN model to the sample alignment and returns the fitted likelihood function"""
    aln = cogent3.load_aligned_seqs(
        "simulated_trees/sim_align/197113_332182_17210.nexus.gz",
        moltype="dna",
    )
    aln = aln.omit_gap_pos(allowed_gap_frac=0)
    tree = cogent3.make_tree("(197113,(332182,17210))")
    model = cogent3.get_app(
        "model", "ssGN", optimise_motif_probs=True, tree=tree, show_progress=False
    )
    result = model(aln)
    return result.lf


def get_seed_tree(num_tips=100):
    # using a midpoint rooted version of the Zhu et al tree
    zhu_tree = cogent3.load_tree(
        "simulated_trees/sim_align/zhu-tree-rerooted-resolved-molclock.nwk"
    )
    # modify tree by extending length of leaf nodes so that the distance from the root to each node
    # is equal (to the initial maximum length).
    tips = zhu_tree.get_tip_names()
    subtips = list(
        rng.choice(list(set(tips).difference({ROOT})), size=num_tips, replace=False)
    ) + [ROOT]
    return zhu_tree.get_sub_tree(subtips)


def sim_alignment(num_seqs=10, align_length=100):
    lf = get_sim_params()
    mprobs = lf.get_motif_probs()
    rates = [p for p in lf.get_param_names() if ">" in p]
    mles = {p: lf.get_param_value(par_name=p) for p in rates}

    sm = cogent3.get_model("ssGN")
    tree = get_seed_tree(num_tips=num_seqs)
    lf = sm.make_likelihood_function(tree)
    lf.set_motif_probs(mprobs)
    for p, v in mles.items():
        lf.set_param_rule(par_name=p, value=v)
    return lf.simulate_alignment(sequence_length=align_length, random_series=rng)


def reroot_at(tip):
    zhu_tree = cogent3.load_tree(
        "simulated_trees/sim_align/zhu-tree-rooted-resolved-molclock.nwk"
    )
    zhu_tree = zhu_tree.unrooted()
    zhu_tree = zhu_tree.rooted_with_tip(tip).bifurcating()
    zhu_tree.write("simulated_trees/sim_align/zhu-tree-rerooted-resolved-molclock.nwk")


if __name__ == "__main__":
    # reroot_at(ROOT)
    taxa = int(sys.argv[1])
    length = int(sys.argv[2])
    sim = sim_alignment(taxa, length)

    with open("simulated_trees/aln_data/" + str(taxa) + ".pkl", "wb") as f:
        pickle.dump(sim, f, pickle.HIGHEST_PROTOCOL)
