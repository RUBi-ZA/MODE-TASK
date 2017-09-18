#!/bin/bash

cd ..

./tsne.py -t tests/pca_test_trj.xtc -p tests/complex.pdb -ag CA -dt rmsd
