#!/usr/bin/env python
"""
Calculates and Returns Diagonals of Correlated Matrix for a given set of modes

meanSquareFluctuation.py
Calculates the mean square fluctuations if carbon atoms:
1) overall msf over all normal modes
2) msf for a specified mode range

Lets say that the user has performed NMA on two coarse grained models of the
same protein and now wants to compare to see if the additional coarse graining
decreased the accuracy. If we obtain the same mean square fluctuations for
each residue then in each model then we can say that the results are comparable
regardless of the coarse graining level. But obviously must compare only the
residues that are common in each model. hence we specify commonResidues (There
is a script that works these out as well so we must just link the output of that)
here i have specified by chain

Author: Caroline Ross: caroross299@gmail.com
August 2017 (Modified by Olivier on November 2020)
"""
import os
import argparse
import numpy as np
from datetime import datetime
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import seaborn as sns
from matplotlib import pyplot as plt
from lib.utils import *

def check_mode_input(text):
    """Checks whether the number of modes is enough to \
    calculate MSF"""
    if not any([len(text.split(":")) < 2 or len(text.split(",")) < 2]):
        print("ERROR: Incorrect number of modes specified")
        sys.exit(1)

def plot(filename, label="", skiprows=1, figsize=(8, 3), linewidth=0.6,
         padding=1, xinterval=3):
    """Generic plotter"""
    dat = pd.read_csv(filename, sep="\t", skiprows=skiprows, index_col=0)
    fig = plt.figure(figsize=figsize, dpi=300)
    ax = plt.subplot(1, 1, 1)
    ax.plot(dat.values, linewidth=linewidth)
    xvals = range(dat.shape[0]) 
    ax.set_xlim(xvals[0]-padding, xvals[-1]+padding)
    ax.scatter(x=xvals, y=dat.values, marker="o", s=8,
               facecolor=(1,1,1,1), edgecolor=(1,0,0), linewidth=linewidth)
    ax.set_xticks(xvals[::xinterval])
    ax.set_xticklabels(dat.index[::xinterval], fontdict={"rotation":90})
    ax.set_xlabel("Residue")
    ax.set_ylabel(label)
    plt.tight_layout()
    plt.savefig("{}.png".format(os.path.splitext(filename)[0]),
                bbox_inches="tight")

def parsePDB(fname, atomT):
    """Parses PDB file"""
    try:	
        f = open(fname, 'r')
        lines = f.readlines()
        f.close()
        f_residues = {}
        for line in lines:
            if line.startswith("ATOM"):
                info = line.split()
                atype = info[2].strip()
                res_type = info[3].strip()
                chain = info[4].strip()
                res = int(info[5].strip())
                if atype == atomT or (atype == "CA" and res_type == "GLY"):
                    if chain in f_residues:
                        if res not in f_residues[chain]:
                            f_residues[chain].append(res)
                    else:
                        f_residues[chain] = [res]
        return f_residues
    except IOError:
        print ("ERROR: File {} not found".format(pdb1))
        sys.exit(1)

def getCommonResidues(c1_residues, c2_residues):
    """Returns common residues shared by 2 structures"""
    common_residues = {}
    for ch in c2_residues:
        for r in c2_residues[ch]:
            if ch in c1_residues:
                if r in c1_residues[ch]:
                    if ch in common_residues:
                        common_residues[ch].append(r)
                    else:
                        common_residues[ch] = [r]
            else:
                break
    return common_residues

def calcMSF(pdb, common_residues, protein_name, wMatrix, vtMatrix,
	    mode_range, specific_residues, atomT):
    """Main function calculating MSF from W and VT matrices"""
    interface_index = []
    CResiduesOrderedByIndex = []
    ResiduesOrderedByIndex = []
    # The pdb file on which NMA analysis was performed, this script handles one
    # model at a time but we can change this
    try:
        f = open(pdb, 'r')
        nma = f.readlines()
        f.close()
    except IOError:
        print ("ERROR: {0} PDB file: {1} not found:".format(protein_name, pdb1))
        sys.exit(1)

    count = 0
    for line in nma:
        if line.startswith("ATOM"):
            info = line.split()
            atype = info[2].strip()
            res_type = info[3].strip()
            chain = info[4].strip()
            res = int(info[5].strip())
            if atype == atomT or (atype == "CA" and res_type == "GLY"):
                if chain in common_residues:
                    if res in common_residues[chain]:
                        # gets the index of the common atoms, as they would
                        # appear in the output W, U and VT matrix. As
                        # these matrices contain info on all atoms
                        CResiduesOrderedByIndex.append(chain+'-'+str(res))
                        interface_index.append(count)

                ResiduesOrderedByIndex.append(chain+'-'+str(res))
                count += 1
    # Specify modes
    total_modes = 0
    try:
        f = open(wMatrix, 'r')
        eigen_values = f.readlines()
        f.close()
    except IOError:
        print ("ERROR: {0} W-matrix file: {1} not found:".format(
               protein_name, wMatrix))
        sys.exit(1)

    total_modes = len(eigen_values)
    res_range = range(count)
    # Create A Full W Inverse Matrix
    # (This is if we want correlation averaged over all modes)
    w_inv = np.zeros((total_modes, total_modes))
    for i in range(total_modes):
        if i > 5:
            w_inv[i, i] = 1 / (float(eigen_values[i].split()[1].strip()))

    print ("INFO: Selecting modes")
    # Create Filtered W Inverse Matrix
    # (This is if we want correlation for a specific mode)
    w_f = np.zeros((total_modes, total_modes))
    for i in mode_range:
        w_f[i, i] = 1 / (float(eigen_values[i].split()[1].strip()))

    # Read In U and VT full Matrix as U is the transpose of VT I only read in
    # VT and create U from the VT matrix.
    print("INFO: Parsing VT_Matrix")
    try:
        fvt = open(vtMatrix, 'r')
        eigen_vectors = fvt.readlines()
        fvt.close()
    except IOError:
        print("ERROR: {0} VT-matrix file: {1} not found".format(
               protein_name, vtMatrix))
        sys.exit(1)

    print("INFO: Transposing VT")
    v_t = np.zeros((total_modes, total_modes))
    u = np.zeros((total_modes, total_modes))
    try:
        for i in range(total_modes):
            vectors = eigen_vectors[i].split()
            for j in range(total_modes):
                vector = float(vectors[j].strip())
                v_t[i, j] = vector
                u[j, i] = vector
        # Calculate Correlation Matrices
        print("INFO: Calculating MSF across all modes")
        w_v_t = np.dot(w_inv, v_t)
        c = np.dot(u, w_v_t)
        print("INFO: Calculating MSF across specific modes")
        # Mode Specific C Matrix
        w_v_tm = np.dot(w_f, v_t)
        CMS = np.dot(u, w_v_tm)
        # Calculate Trace of the Correlation Matrices
        trace_c = np.zeros((total_modes//3, total_modes//3))
        trace_c_m = np.zeros((total_modes//3, total_modes//3))
        for i in range(0, total_modes, 3):
            for j in range(0, total_modes, 3):
                trace = 0
                trace_m = 0
                for k in range(3):
                    trace = trace + c[i + k, j + k]
                    trace_m = trace_m + CMS[i + k, j + k]
                trace_c[i//3, j//3] = trace
                trace_c_m[i//3, j//3] = trace_m
        # Print the diagonal values per residue
        
        # Process/plot RMSF from all modes
        outfilename_allmodes = args.outdir + "/" + protein_name+ "_msf.txt"
        w = open(outfilename_allmodes, 'w')
        w.write("MSF Calculated for "+pdb)
        w.write("\nRes\tm.s.f\n")
        for i in res_range:
            w.write(str(ResiduesOrderedByIndex[i]) + "\t" + str(trace_c[i, i]) + "\n")
        w.close()
        plot(outfilename_allmodes, label="MSF from all modes", skiprows=1)
        
        # Process/plot RMSF from specified modes
        outfilename_selmodes = args.outdir + "/" + protein_name +"_msfSpecificModes.txt" 
        w = open(outfilename_selmodes, 'w')
        w.write("MSF Calculated for "+pdb)
        w.write("\n**********************\nMODES SELECTED FOR CALCULATIONS:\n")
        for m in mode_range:
            w.write(str(m + 1) + " ; ")
        w.write("\n**********************")
        w.write("\nRes\tm.s.f\n")
        for i in res_range:
            w.write(str(ResiduesOrderedByIndex[i])+ "\t" + str(trace_c_m[i, i]) + "\n")
        w.close()
        plot(outfilename_selmodes, label="MSF from  modes {}".format(args.modes), skiprows=5)

        if specific_residues:
            w = open(args.outdir + "/" + protein_name + "CommonResidues_msf.txt", 'w')
            w.write("MSF Calculated for "+pdb)
            w.write("\nRes\tm.s.f\n")
            for k,i in enumerate(interface_index):
                w.write(str(CResiduesOrderedByIndex[k]) + "\t" + str(trace_c[i, i]) + "\n")
            w.close()
            w = open(args.outdir + "/" + protein_name +"_CommonResidues_msfSpecificModes.txt", 'w')
            w.write("MSF Calculated for "+pdb)
            w.write("\n**********************\nMODES SELECTED FOR CALCULATIONS:\n")
            for m in mode_range:
                w.write(str(m+1)+" ; ")
            w.write("\n**********************")
            w.write("\nRes\tm.s.f\n")
            for k, i in enumerate(interface_index):
                w.write(str(CResiduesOrderedByIndex[k]) + "\t" + str(trace_c_m[i, i]) + "\n")
            w.close()
    except IndexError:
        print ("ERROR: PDB file and ANM matrices are incompatible. \
                Check for 1. incorrect pdb file, matrices or specified modes")
        sys.exit(1)

def main(args):
    """Main program"""
    atomT = args.atomType.upper()
    if atomT != 'CA' and atomT != 'CB':
        print ("ERROR: Unrecognised atom type (Option: CA or CB")
        sys.exit(1)
    specificResidues = True
    pdb1 = args.pdb
    if args.pdbC == "none":
        specificResidues = False
    else:
        pdb2 = args.pdbC
    modes = args.modes
    check_mode_input(modes)
    mode_range = []
    if ':' in modes:
        mRange = modes.split(':')
        mode_range = range(int(mRange[0])-1,int(mRange[1]))
    else:
        modesSpecific = modes.split(',')
        for m in modesSpecific:
            mode_range.append(int(m)-1)
    if specificResidues:
        pdb1_residues = parsePDB(pdb1,atomT)
        pdb2_residues = parsePDB(pdb2,atomT)
        common_residues = getCommonResidues(pdb1_residues,pdb2_residues)
	# Process 1st protein
        protein_name1 = "PDB1"
        wMatrix1 = args.wMatrix
        vtMatrix1 = args.vtMatrix
        calcMSF(pdb1, common_residues, protein_name1, wMatrix1, vtMatrix1,
                mode_range, specificResidues, atomT)
	# Process 2nd protein
        protein_name2 = "PDBCompare"
        wMatrixC = args.wMatrixC
        vtMatrixC = args.vtMatrixC
        calcMSF(pdb2, common_residues, protein_name2, wMatrixC, vtMatrixC,
                mode_range, specificResidues, atomT)
        #write a list of common residues to text file
        w = open(args.outdir + "/common_residues_msf", 'w')
        w.write("Residues Common in Both PDB Files: "+pdb1+" and "+pdb2+"\n")
        w.write(str(common_residues))
        w.close()
    else:
        protein_name1 = "PDB1"
        pdb1_residues = parsePDB(pdb1,atomT)
        wMatrix1 = args.wMatrix
        vtMatrix1 = args.vtMatrix
        calcMSF(pdb1, pdb1_residues, protein_name1, wMatrix1, vtMatrix1,
		mode_range, specificResidues, atomT)

silent = False
stream = sys.stdout

def log(message):
    global silent
    global stream
    if not silent:
        print_err(message)

if __name__ == "__main__":
    sns.set_style("whitegrid")
    plt.rcParams["font.family"] = "serif"
    plt.rcParams['xtick.labelsize'] = 6
    plt.rcParams['ytick.labelsize'] = 6
    parser = argparse.ArgumentParser()
    # standard arguments for logging
    parser.add_argument("--silent", action='store_true',
			default=False, help="Turn off logging")
    parser.add_argument("--welcome", default="true",
			help="Display welcome message (true/false)")
    parser.add_argument("--log-file", default=None,
			help="Output log file (default: standard output)")
    parser.add_argument("--outdir", help="Output directory", default="output")
    # custom arguments
    parser.add_argument("--pdb", help="Input")
    parser.add_argument("--pdbC", default="none",
			help="Enter a second PDB file to compare MSF between two models")
    parser.add_argument("--modes", default="7:27",
			help="Enter a select range of modes in format M1:M2 \
			E.g To calculate MSF over the first 20 non-zero modes \
			enter --modes 7:27 OR Calculate the MSF for a \
			combination of specific modes Enter mode numbers \
			separated by a comma Eg: --modes 1,5,7")
    parser.add_argument("--wMatrix",
			help="Text file of Eigevalues of pdb, in format \
			output from ANM.cpp")
    parser.add_argument("--vtMatrix",
			help="Text file of Eigevectors of pdb, in row (VT)\
			 format output from ANM.cpp")
    parser.add_argument("--wMatrixC", default="none",
			help="Text file of Eigevalues of pdbC, in format \
			output from ANM.cpp")
    parser.add_argument("--vtMatrixC", default="none",
			help="Text file of Eigevectors of pdbC, in row (VT) \
			format output from ANM.cpp")
    parser.add_argument("--atomType", default='X',
			help="Enter CA to select alpha carbons or CB to select \
			 beta carbons")
    args = parser.parse_args()

    if args.welcome == "true":
        welcome_msg("Mean square fluctuations",
		    "Caroline Ross (caroross299@gmail.com)")
    # Check if args supplied by user
    if len(sys.argv) > 1:
        if not os.path.isdir(args.outdir):
            os.makedirs(args.outdir)
        # set up logging
        silent = args.silent
        if args.log_file:
            stream = open(args.log_file, 'w')
        start = datetime.now()
        log("Started at: %s" % str(start))
        main(args)
        end = datetime.now()
        time_taken = format_seconds((end - start).seconds)
        log("Completed at: %s" % str(end))
        log("- Total time: %s" % str(time_taken))
        # close logging stream
        stream.close()
    else:
        print("ERROR: No arguments provided. Use -h to view help")
        sys.exit(1)
