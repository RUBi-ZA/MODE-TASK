#!/usr/bin/env python
"""Wrapper script/plotter for ANM calculation"""
__author__ = ["Michael Glenister", "Olivier Sheik Amamuddy"]
__copyright__ = "Copyright 2020, Research Unit in Bioinformatics"
__license__ = "GNU GPL 3.0"
__version__ = "1.1"
__email__ = "oliserand@gmail.com"
__status__ = "Production"
__date__ = "25th Nov 2020"

import os
import sys
import argparse
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import seaborn as sns
from matplotlib import pyplot as plt

def plot(filename, nmodes=50, figsize=(8, 3)):
    """Plots the 1 n eigenvalues"""
    dat = pd.read_csv(filename, sep=r"\s+", names=["res", "eigval"], header=None)
    dat = dat[:nmodes]
    yeps = (dat.eigval.max() - dat.eigval.min())/200
    fig = plt.figure(figsize=figsize, dpi=300)
    axis = sns.barplot(x="res", y="eigval", data=dat, color="red")
    _ = [axis.text(x=i, y=dat.eigval.values[i]+yeps,
		   s=str(dat.eigval[i].round(6)),
                   fontdict={"fontsize":7, "rotation":90, "ha": "center",
                             "va": "bottom"})
	 for i in dat.index]
    axis.set_xticklabels(dat.res,
                       fontdict={"rotation":90})
    axis.set_ylim(0, dat.eigval.max()*1.15)
    axis.set_xlabel("Mode number",
                    fontdict={"fontweight":"semibold", "fontsize":8})
    axis.set_ylabel("Eigenvalue",
                    fontdict={"fontweight":"semibold", "fontsize":8})
    axis.set_title("First {} eigenvalues".format(nmodes),
                   fontdict={"fontweight":"bold"})
    plt.tight_layout()
    sns.despine(fig)
    outfilename = "{}".format(os.path.join(os.path.dirname(filename),
                                  "eigenvectors.png"))
    plt.savefig(outfilename, bbox_inches="tight")
    print("INFO: Wrote {}".format(outfilename))

def parse_args():
    """Argument parser"""
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb", help="PDB input file")
    parser.add_argument("--cutoff", default="15",
                        help="Cutoff radius in Angstroms (default: 15)")
    parser.add_argument("--atomType", help="Select 'CA' or 'CB' atoms")
    parser.add_argument("--outdir", default="output",
                        help="Output directory")
    args = parser.parse_args()
    return args

def main(args):
    """Program main"""
    # Check if required directories exist
    if not os.path.isdir(args.outdir):
        os.makedirs(args.outdir)
    if len(sys.argv) > 1:
        if not os.path.isfile("ANM"):
            os.system("g++ -I cpp/src/ ANM.cpp -o ANM")
            if not os.path.isfile("ANM"):
                print ("ERROR: 'ANM' executable not found.")
                print ("Please ensure it exists in the same directory as this script.")
                sys.exit(1)
        # Run getEigenVectors executable
        os.system("./ANM --pdb {0} --cutoff {1} --atomType {2} --outdir {3}".format(
                  args.pdb, args.cutoff, args.atomType, args.outdir))
        print("INFO: Wrote ANM files in {}".format(args.outdir))
        plt.rcParams["font.size"] = 7
        plt.rcParams["font.family"] = "serif"
        outfilename = "{}".format(os.path.join(args.outdir, "W_values.txt"))
        plot(filename=outfilename)
    else:
        print('ERROR: No arguments provided. Use -h to view help')

if __name__ == "__main__":
    main(parse_args())
