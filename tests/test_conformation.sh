#!/bin/bash

file="input_pdb"
while IFS= read line
do
	prot=${line}
done <${file}

cd ..

./conformationMode.py --pdbConf output/aligned.pdb --pdbANM output/ComplexCG.pdb --vtMatrix output/VT_values.txt
