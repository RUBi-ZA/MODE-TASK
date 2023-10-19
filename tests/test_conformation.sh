#!/bin/bash

file="input_pdb"
while IFS= read line
do
	prot=${line}
done <${file}

cd ../src

./conformationMode.py --pdbConf ../output/ComplexCG4.pdb --pdbANM ../output/ComplexCG4.pdb --vtMatrix ../output/VT_values.txt --atomType CB
