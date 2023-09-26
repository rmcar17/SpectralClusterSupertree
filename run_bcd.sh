#!/bin/sh
time {
java -jar bcd/BCDSupertrees.jar -L OFF -s true -f NEWICK -d NEWICK -o tmp.nwk -t 1 -B "$1" &> /dev/null
wait
cat tmp.nwk
wait
rm tmp.nwk
}