#!/bin/bash
cd ..

g++ -I cpp/input/ getEigenVectors.cpp -o getEigenVectors

./getEigenVectors --vt output/VT_values.txt --mode 1
