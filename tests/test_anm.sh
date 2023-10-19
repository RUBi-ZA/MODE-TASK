#!/bin/bash

file="input_pdb"
while IFS= read line
do
    prot=${line}
done <${file}

cd ../src

./coarseGrain.py --pdb ../tests/${prot}.pdb --cg 4 --atomType CB --outdir ../output/

g++ -I cpp/src/ ANM.cpp -o ANM

../src/ANM --pdb ../output/ComplexCG4.pdb --cutoff 24 --atomType CB







