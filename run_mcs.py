import sys
from min_cut_supertree import min_cut_supertree
from run_scs import parse_trees

if __name__ == "__main__":
    input_trees = parse_trees(sys.argv[1])
    supertree = min_cut_supertree(input_trees)
    print(supertree)
