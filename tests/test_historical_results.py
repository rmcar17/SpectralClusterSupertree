from helpers import load_expected_tree_file, load_source_tree_file, scs_test


def test_dcm_iq():
    expected = load_expected_tree_file("dcm_iq_expected.tre")
    source_trees = load_source_tree_file("dcm_iq_source.tre")

    scs_test(source_trees, expected, pcg_weighting="branch")


def test_smidgen_og():
    expected = load_expected_tree_file("smidgen_og_expected.tre")
    source_trees = load_source_tree_file("smidgen_og_source.tre")

    scs_test(source_trees, expected, pcg_weighting="branch")


def test_smidgen_og_5500():
    expected = load_expected_tree_file("smidgen_og_5500_expected.tre")
    source_trees = load_source_tree_file("smidgen_og_5500_source.tre")

    scs_test(source_trees, expected, pcg_weighting="branch")


def test_supertriplets():
    expected = load_expected_tree_file("supertriplets_expected.tre")
    source_trees = load_source_tree_file("supertriplets_source.tre")

    scs_test(source_trees, expected, pcg_weighting="depth")
