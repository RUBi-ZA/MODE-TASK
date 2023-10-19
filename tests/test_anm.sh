#!/bin/bash

cd ../src

g++ -I cpp/src/ ANM.cpp -o ANM

../src/ANM --pdb ../outdir/ComplexCG4.pdb --cutoff 24 --atomType CB
