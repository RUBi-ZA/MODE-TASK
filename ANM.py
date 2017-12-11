#!/usr/bin/env python

# Simple python script to call getEigenVectors executable

import os
import sys
import argparse

# parse cmd arguments
parser = argparse.ArgumentParser()

# custom arguments
parser.add_argument("--pdb", help="PDB input file")
parser.add_argument("--cutoff", help="Cuttoff radius in Angstroms. Default: 15", default="15")
parser.add_argument("--atomType", help="Select 'CA' or 'CB' atoms")
parser.add_argument("--outdir", help="Directory to generate output to", default="output")

args = parser.parse_args()

# Check if required directories exist
if not os.path.isdir(args.outdir):
    os.makedirs(args.outdir)

# Check if args supplied by user
if len(sys.argv) > 1:

    # Check if getEigenVectors exists
    if not os.path.isfile("ANM"):
        os.system("g++ -I cpp/src/ ANM.cpp -o ANM")
        if not os.path.isfile("ANM"):
            print ("Error: 'ANM' executable not found.")
            print ("Please ensure it exists in the same directory as this script.")
            sys.exit()

    # Run getEigenVectors executable
    os.system("./ANM --pdb " + args.pdb + " --cutoff " + args.cutoff + " --atomType " + args.atomType + " --outdir " + args.outdir)

else:
    print ('No arguments provided. Use -h to view help')
