#!/bin/bash

cd ../src

g++ -I cpp/src/ ANM.cpp -o ANM

./ANM --pdb output/ComplexCG4.pdb --cutoff 24 --atomType CB
