#!/bin/bash

file="input_pdb"
while IFS= read line
do
	prot=${line}
done <${file}

cd ../src

./coarseGrain.py --pdb tests/${prot}.pdb --cg 4 --atomType CB
