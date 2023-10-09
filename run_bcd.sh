#!/bin/sh
time {
java -jar bcd/BCDSupertrees.jar -L OFF -s true -w BRANCH_LENGTH -f NEWICK -d NEWICK -o tmp.nwk -t 1 -B "$1" &> /dev/null
wait # weighting (-w) can either be UNIT_WEIGHT, TREE_WEIGHT, BRANCH_LENGTH or mroe
cat tmp.nwk
wait
rm tmp.nwk
}