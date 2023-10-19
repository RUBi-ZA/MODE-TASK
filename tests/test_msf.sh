#!/bin/bash

cd ../src

./meanSquareFluctuations.py --pdb output/ComplexCG4.pdb --wMatrix output/W_values.txt --vtMatrix output/VT_values.txt --atomType CB
