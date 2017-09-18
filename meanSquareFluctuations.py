#!/usr/bin/env python

# meanSquareFluctuation.py
# Calcualtes the mean square fluctutaions if carbon atoms:
# 1) overall msf over all normal modes
# 2) msf for a specified mode range
# Author: Caroline Ross: caroross299@gmail.com
# August 2017

# Calculates and Returns Diagonals of Correlated Matrix for a given set of modes

import os
import argparse
from datetime import datetime

from lib.utils import *
import numpy as np


# Lets say that the user has performed NMA on two coarse grained models of the same protein and now wants to compare
# to see if the additional coarse graining decreased the accuracy. If we obtain the same mean square fluctuations for
# each residue then in each model then we can say that the results are comparable regardless of the coarse graining
# level. But obviously must compare only the residues that are common in each model. hence we specify commonResidues
# (There is a script that works these out as well so we must just link the output of that) here i have specified by
# chain

def main(args):

    atomT = args.atomType.upper()
    if atomT!='CA' and atomT!='CB':
	print '\n**************************************\nUnrecognised atom type\nInput Options:\nCA: to select alpha carbon atoms\nCB: to select beta carbon atoms\n**************************************'
	sys.exit()

    # residues common to multiple pdb files (if user has input two conformations)
    specificResidues = True
    pdb1 = args.pdb
    if args.pdbConf2 == "none":
	specificResidues = False
	pdb2 = args.pdb
    else:
	pdb2=args.pdbConf2

    #getCommonResidues:
    # Takes two pdb models and determines the common residues
    #####################################
    try:	
        f = open(pdb1, 'r')
        lines_c1 = f.readlines()
        f.close()
    except IOError:
        print '\n**************************************\nFILE '+pdb1+' NOT FOUND:\n**************************************\n'
	sys.exit()
    try:
        f = open(pdb2, 'r')
        lines_c2 = f.readlines()
        f.close()
    except IOError:
        print '\n**************************************\nFILE '+pdb2+' NOT FOUND:\n**************************************\n'
	sys.exit()
    #####################################

    c1_residues = {}
    c2_residues = {}

    common_residues = {}

    for line in lines_c1:
        if line.startswith("ATOM"):
            info = line.split()
            atype = info[2].strip()
            res_type = info[3].strip()
            chain = info[4].strip()
            res = int(info[5].strip())
            if atype == atomT or (atype == "CA" and res_type == "GLY"):
                if chain in c1_residues:
                    if res not in c1_residues[chain]:
                        c1_residues[chain].append(res)
                else:
	            c1_residues[chain] = [res]

    for line in lines_c2:
        if line.startswith("ATOM"):
            info = line.split()
            atype = info[2].strip()
            res_type = info[3].strip()
            chain = info[4].strip()
            res = int(info[5].strip())
            if atype == atomT or (atype == "CA" and res_type == "GLY"):
                if chain in c2_residues:
                    if res not in c2_residues[chain]:
                        c2_residues[chain].append(res)
                else:
	            c2_residues[chain] = [res]

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
    
    w = open(args.outdir + "/common_residues_msf", 'w')
    w.write(str(common_residues))
    w.close()

    interface_index = []
    CResiduesOrderedByIndex = []
    ResiduesOrderedByIndex = []

    # The pdb file on which NMA analysis was performed, this script handles one model at a time but we can change this
    try:
        f = open(pdb1, 'r')
        nma = f.readlines()
        f.close()
    except IOError:
        print '\n**************************************\nFILE '+pdb1+' NOT FOUND:\n**************************************\n'
	sys.exit()

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
                        # gets the index of the common atoms, as they would appear in the output W, U and VT matrix. As
                        # these matrices contain info on all atoms
		        CResiduesOrderedByIndex.append(chain+'-'+str(res))
                        interface_index.append(count)
		
		ResiduesOrderedByIndex.append(chain+'-'+str(res))
                count += 1
  
    protein_name = args.pdb
    protein_name = protein_name[protein_name.rfind("/")+1:protein_name.rfind("/")+5]

    # Specify modes
    total_modes = 0
    try: 
        f = open(args.wMatrix, 'r')
        eigen_values = f.readlines()
        f.close()
    except IOError:
        print '\n**************************************\nFILE '+args.wMatrix+' NOT FOUND:\n**************************************\n'
	sys.exit()
    for i in range(len(eigen_values)):
        total_modes += 1

    first_mode = args.firstMode  
    last_mode = args.lastMode  

    # If user fails to provide first and last modes get default values
    if args.firstMode == 0 and args.lastMode == 0:
	last_mode = total_modes - 6
        first_mode = last_mode - 19
        

    # Llama
    first_res = 0
    # Specify Residue Indexes
    # Get first residue number
    '''for line in nma:  # Future work to include this in existing pdb file read above
        if line.startswith("ATOM"):
            if first_res == 0:
                info = line.split()
                first_res = int(info[1].strip())'''

    res_range = range(count)
    print "Residue count: " + str(count)

    mode_range = range(first_mode-1, last_mode)
    

    # Create A Full W Inverse Matrix (This is if we want correlation averaged over all modes)
    w_inv = np.zeros((total_modes, total_modes))
    for i in range(total_modes):
        if i < total_modes - 6:
            w_inv[i, i] = 1 / (float(eigen_values[i].split()[1].strip()))
    
    print "Selecting modes"
    # Create Filtered W Inverse Matrix (This is if we want correlation for a specific mode)
    w_f = np.zeros((total_modes, total_modes))
    for i in mode_range:
        w_f[i, i] = 1 / (float(eigen_values[i].split()[1].strip()))

    # Read In U and VT full Matrix as U is the transpose of VT I only read in VT and create U from the VT matrix
    # info. So we can exclude U output from C++ script for faster analysis
    print "Parsing VT_Matrix"

    try:
        fvt = open(args.vtMatrix, 'r')  
        eigen_vectors = fvt.readlines()
        fvt.close()
    except IOError:
        print '\n**************************************\nFILE '+args.vtMatrix+' NOT FOUND:\n**************************************\n'
	sys.exit()

    print "Calculating Transpose of VT"
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
	print "Calculating MSF across all modes"
	w_v_t = np.dot(w_inv, v_t)
	# print "Correlations Calculated"
	c = np.dot(u, w_v_t)
	# print "Correlations Calculated"
    
	print "Calculating MSF across modes: "+str(first_mode)+"-"+str(last_mode)	
	# Mode Specific C Matrix
	w_v_tm = np.dot(w_f, v_t)
	CMS = np.dot(u, w_v_tm)
    
    
	# Calculate Trace of the Correlation Matrices
	trace_c = np.zeros((total_modes / 3, total_modes / 3))
	trace_c_m = np.zeros((total_modes / 3, total_modes / 3))
    
	for i in range(0, total_modes, 3):
	    for j in range(0, total_modes, 3):
		trace = 0
		trace_m = 0
		for k in range(3):
		    trace = trace + c[i + k, j + k]
		    trace_m = trace_m + CMS[i + k, j + k]
		trace_c[i / 3, j / 3] = trace
		trace_c_m[i / 3, j / 3] = trace_m
    
	# Print the diagonal values per residue
	w = open(args.outdir + "/" + protein_name+ "_msf.txt", 'w')
	w.write("Res\tm.s.f\n")
     
	for i in res_range:
	    w.write(str(ResiduesOrderedByIndex[i]) + "\t" + str(trace_c[i, i]) + "\n")
	w.close()
    
	w = open(args.outdir + "/" + protein_name +"_msfModes"+ str(first_mode)+"_"+str(last_mode)+".txt", 'w')
	w.write("Res\tm.s.f\n")
	for i in res_range:
	    w.write(str(ResiduesOrderedByIndex[i])+ "\t" + str(trace_c_m[i, i]) + "\n")
	w.close()
	
      
	if specificResidues:
	    w = open(args.outdir + "/" + protein_name + "CommonResidues_msf.txt", 'w')
	    w.write("Res\tm.s.f\n")
	    for k,i in enumerate(interface_index): 
	      
		w.write(str(CResiduesOrderedByIndex[k]) + "\t" + str(trace_c[i, i]) + "\n")
	    w.close()
    
	    w = open(args.outdir + "/" + protein_name +"_CommonResidues_msfModes"+ str(first_mode)+"_"+str(last_mode)+".txt", 'w')	
	    w.write("Res\tm.s.f\n")
	    for k,i in enumerate(interface_index):
		w.write(str(CResiduesOrderedByIndex[k]) + "\t" + str(trace_c_m[i, i]) + "\n")
	    w.close()
    except IndexError:
	print '\n**************************************\nERROR!!\nINPUT FILES ARE NOT COMPATIBLE\n**************************************\n'
        sys.exit()	
	

silent = False
stream = sys.stdout


def log(message):
    global silent
    global stream

    if not silent:
        print >> stream, message


if __name__ == "__main__":
    # parse cmd arguments
    parser = argparse.ArgumentParser()

    # standard arguments for logging
    parser.add_argument("--silent", help="Turn off logging", action='store_true', default=False)
    parser.add_argument("--welcome", help="Display welcome message (true/false)", default="true")
    parser.add_argument("--log-file", help="Output log file (default: standard output)", default=None)
    parser.add_argument("--outdir", help="Output directory", default="output")

    # custom arguments
    parser.add_argument("--pdb", help="Input") 
    parser.add_argument("--pdbConf2", help="", default = "none")  
    parser.add_argument("--firstMode", help="[int]", default=0, type=int) #default = total modes-25
    parser.add_argument("--lastMode", help="[int]", default=0, type=int)  # use last non-trivial mode if no user input provided
    parser.add_argument("--wMatrix", help="W matrix input file that was output from C++ Scripts")
    parser.add_argument("--vtMatrix", help="U and VT full Matrix")
    parser.add_argument("--atomType", help="Enter CA to select alpha carbons or CB to select beta carbons", default='CA')

    args = parser.parse_args()

    if args.welcome == "true":
        welcome_msg("Mean square fluctuations", "Caroline Ross (caroross299@gmail.com)")

    # Check if args supplied by user
    if len(sys.argv) > 1:
        # Check modes
        if args.firstMode > args.lastMode:
            print '\n**************************************\nFirst mode cannot be greater than last mode\n**************************************\n'
            sys.exit()
        if args.firstMode<0 or args.lastMode<0:
	    print '\n**************************************\nINVALID INPUT\nFirst and last mode must be positive integers\n**************************************\n'
            sys.exit()


        # Check if required directories exist
        if not os.path.isdir(args.outdir):
            os.makedirs(args.outdir)

        # set up logging
        silent = args.silent

        if args.log_file:
            stream = open(args.log_file, 'w')

        start = datetime.now()
        log("Started at: %s" % str(start))

        # run script
        main(args)

        end = datetime.now()
        time_taken = format_seconds((end - start).seconds)

        log("Completed at: %s" % str(end))
        log("- Total time: %s" % str(time_taken))

        # close logging stream
        stream.close()
    else:
        print "No arguments provided. Use -h to view help"
