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


def parsePDB(fname,atomT):
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
        print ('\n**************************************\nFILE ' + pdb1 + ' NOT FOUND:\n**************************************\n')
        sys.exit()
        
        
def getCommonResidues(c1_residues,c2_residues):
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


def calcMSF(pdb,common_residues,protein_name,wMatrix,vtMatrix,mode_range,sr,atomT):
    
    interface_index = []
    CResiduesOrderedByIndex = []
    ResiduesOrderedByIndex = []

    # The pdb file on which NMA analysis was performed, this script handles one model at a time but we can change this
    try:
        f = open(pdb, 'r')
        nma = f.readlines()
        f.close()
    except IOError:
        print ('\n**************************************\n'+protein_name+' PDB FILE: '+pdb1+' NOT FOUND:\n**************************************\n')
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
  
    # Specify modes
    total_modes = 0
    try: 
        f = open(wMatrix, 'r')
        eigen_values = f.readlines()
        f.close()
    except IOError:
        print ('\n**************************************\n'+protein_name+' W-MATRIX FILE: '+wMatrix+' NOT FOUND:\n**************************************\n')
        sys.exit()
	
    total_modes = len(eigen_values)
    res_range = range(count)
    #Create A Full W Inverse Matrix (This is if we want correlation averaged over all modes)
    w_inv = np.zeros((total_modes, total_modes))
    for i in range(total_modes):
        if i >5:
            w_inv[i, i] = 1 / (float(eigen_values[i].split()[1].strip()))

    print ("Selecting modes")
    # Create Filtered W Inverse Matrix (This is if we want correlation for a specific mode)
    w_f = np.zeros((total_modes, total_modes))
    for i in mode_range:
        w_f[i, i] = 1 / (float(eigen_values[i].split()[1].strip()))


    # Read In U and VT full Matrix as U is the transpose of VT I only read in VT and create U from the VT matrix
    # info. So we can exclude U output from C++ script for faster analysis
    print ("Parsing VT_Matrix")

    try:
        fvt = open(vtMatrix, 'r')  
        eigen_vectors = fvt.readlines()
        fvt.close()
    except IOError:
        print ('\n**************************************\n'+protein_name+' VT-MATRIX FILE: '+vtMatrix+' NOT FOUND:\n**************************************\n')
        sys.exit()

    print ("Calculating Transpose of VT")
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
        print ("Calculating MSF across all modes")
        w_v_t = np.dot(w_inv, v_t)
        # print "Correlations Calculated"
        c = np.dot(u, w_v_t)
        # print "Correlations Calculated"
    
        print ("Calculating MSF across modes specific modes")	
        # Mode Specific C Matrix
        w_v_tm = np.dot(w_f, v_t)
        CMS = np.dot(u, w_v_tm)
    
    
        # Calculate Trace of the Correlation Matrices
        trace_c = np.zeros((total_modes // 3, total_modes // 3))
        trace_c_m = np.zeros((total_modes // 3, total_modes // 3))
    
        for i in range(0, total_modes, 3):
            for j in range(0, total_modes, 3):
                trace = 0
                trace_m = 0
                for k in range(3):
                    trace = trace + c[i + k, j + k]
                    trace_m = trace_m + CMS[i + k, j + k]
                trace_c[i // 3, j // 3] = trace
                trace_c_m[i // 3, j // 3] = trace_m
    	    
    	    
        # Print the diagonal values per residue
        w = open(args.outdir + "/" + protein_name+ "_msf.txt", 'w')
        w.write("MSF Calculated for "+pdb)
        w.write("\nRes\tm.s.f\n")
     
        for i in res_range:
            w.write(str(ResiduesOrderedByIndex[i]) + "\t" + str(trace_c[i, i]) + "\n")
        w.close()
        
        
        w = open(args.outdir + "/" + protein_name +"_msfSpecificModes.txt", 'w')
        w.write("MSF Calculated for "+pdb)
        w.write("\n**********************\nMODES SELECTED FOR CALCULATIONS:\n")
        for m in mode_range:
            w.write(str(m+1)+" ; ")
        w.write("\n**********************")
        w.write("\nRes\tm.s.f\n")
        for i in res_range:
            w.write(str(ResiduesOrderedByIndex[i])+ "\t" + str(trace_c_m[i, i]) + "\n")
        w.close()    
        
        
      
        if sr:
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
            for k,i in enumerate(interface_index):
                w.write(str(CResiduesOrderedByIndex[k]) + "\t" + str(trace_c_m[i, i]) + "\n")
            w.close()

    except IndexError:
        print ('\n**************************************\nERORR!! PDB FILE AND ANM MATRICES ARE IMCOMPATABLE\nCHECK INPUT PARAMETERS FOR:\n1) INCORRECT PDB FILE \n2) INCORRECT MATRICES \n3) INCORRECT SPECIFIED MODES\n**************************************\n')
        sys.exit()


def main(args):

    atomT = args.atomType.upper()
    if atomT!='CA' and atomT!='CB':
        print ('\n**************************************\nUnrecognised atom type\nInput Options:\nCA: to select alpha carbon atoms\nCB: to select beta carbon atoms\n**************************************')
        sys.exit()

    specificResidues = True
    pdb1 = args.pdb
    if args.pdbC == "none":
        specificResidues = False
    else:
        pdb2=args.pdbC

    modes = args.modes	
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
        protein_name1 = "PDB1"
        wMatrix1 = args.wMatrix
        vtMatrix1 = args.vtMatrix 
        calcMSF(pdb1,common_residues,protein_name1,wMatrix1,vtMatrix1,mode_range,specificResidues,atomT)
	
	
        protein_name2 = "PDBCompare"
        wMatrixC = args.wMatrixC
        vtMatrixC = args.vtMatrixC 
        calcMSF(pdb2,common_residues,protein_name2,wMatrixC,vtMatrixC,mode_range,specificResidues,atomT)
	
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
        calcMSF(pdb1,pdb1_residues,protein_name1,wMatrix1,vtMatrix1,mode_range,specificResidues,atomT)
        
        
        




	

silent = False
stream = sys.stdout


def log(message):
    global silent
    global stream

    if not silent:
        print_err(message)


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
    parser.add_argument("--pdbC", help="Enter a second PDB file to compare MSF bewteen two models", default = "none")  
    parser.add_argument("--modes", help="Enter a select range of modes in format M1:M2\nE.g To calculate MSF over the first 20 non-zero modes enter --modes 7:27\nOR\nCalculate the MSF for a combination of specific modes\nEnter mode numbers separated by a comma\nEg: --modes 1,5,7", default="7:27") 
    parser.add_argument("--wMatrix", help="Text file of Eigevalues of pdb, in format output from ANM.cpp")
    parser.add_argument("--vtMatrix", help="Text file of Eigevectors of pdb, in row (VT) format output from ANM.cpp")
    parser.add_argument("--wMatrixC", help="Text file of Eigevalues of pdbC, in format output from ANM.cpp", default = "none")
    parser.add_argument("--vtMatrixC", help="Text file of Eigevectors of pdbC, in row (VT) format output from ANM.cpp", default = "none")
    parser.add_argument("--atomType", help="Enter CA to select alpha carbons or CB to select beta carbons", default='CA')

    args = parser.parse_args()

    if args.welcome == "true":
        welcome_msg("Mean square fluctuations", "Caroline Ross (caroross299@gmail.com)")

    # Check if args supplied by user
    if len(sys.argv) > 1:

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
        print ('No arguments provided. Use -h to view help')
