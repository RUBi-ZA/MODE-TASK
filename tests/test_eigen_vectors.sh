#!/bin/bash
cd ..

g++ -I cpp/input/ getEigenVectors.cpp -o getEigenVectors

./getEigenVectors --vtMatrix output/VT_values.txt --mode 1
