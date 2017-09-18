#!/bin/bash

cd ..

g++ -I cpp/src/ ANM.cpp -o ANM

./ANM --pdb output/ComplexCG.pdb --cutoff 24
