#!/usr/bin/env python

# Simple python script to call getEigenVectors executable

import os
import sys
import argparse

# parse cmd arguments
parser = argparse.ArgumentParser()

# custom arguments
parser.add_argument("--vt", help="VT matrix file")
parser.add_argument("--mode", help="Specify the index of the mode you wish to target")
parser.add_argument("--direction", help="Direction of overlap correction (Default -1)", default="-1")
parser.add_argument("--outdir", help="Directory to generate output to", default="output")

args = parser.parse_args()

# Check if required directories exist
if not os.path.isdir(args.outdir):
    os.makedirs(args.outdir)

# Check if args supplied by user
if len(sys.argv) > 1:

    # Check if getEigenVectors exists
    if not os.path.isfile("getEigenVectors"):
        os.system("g++ -I cpp/input/ getEigenVectors.cpp -o getEigenVectors")
        if not os.path.isfile("getEigenVectors"):
            print ("Error: 'getEigenVectors' executable not found.")
            print ("Please ensure it exists in the same directory as this script.")
            sys.exit()

    # Run getEigenVectors executable
    os.system("./getEigenVectors --vt " + args.vt + " --mode " + args.mode + " --direction " + args.direction + " --outdir " + args.outdir)

else:
    print ('No arguments provided. Use -h to view help')
