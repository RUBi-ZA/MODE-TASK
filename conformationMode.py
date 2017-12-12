#!/usr/bin/env python

# coformationalMode.py
# Identifies normal modes that act in the direction of a conformational change
# Author: Caroline Ross: caroross299@gmail.com
# August 2017


import os
import sys

# Append to system path to allow for importing local libraries
pathname = os.path.dirname(sys.argv[0])
sys.path.append(os.path.abspath(pathname))

import argparse
from datetime import datetime
from lib.utils import *
import numpy as np
from math import sqrt
import numpy as np

from lib import sdrms

def calcConformation(delta_r,vtMatrix,nma_index,outdir,outfile):
    # Calculate the magnitude
    correlationDR = []
    csum=0
    mag_d_r = 0
    for i,r in enumerate(delta_r):
        mag_d_r += r * r
        csum = csum+r*r
        if i%3==2:
            csum = sqrt(csum)
            correlationDR.append(csum)
            csum=0
    mag_d_r = sqrt(mag_d_r)


    try:   
        f = open(vtMatrix, 'r')
        vectors = f.readlines()
        f.close()
    except IOError:
        print ('\n**************************************\nFILE '+vtMatrix+' NOT FOUND:\n**************************************\n')
        sys.exit()
    mode_range = range(len(vectors)-6) #excludes the trivial modes
 
    
    output = {}
    try:
        for mode in mode_range:
            correlationMode = []
            overlap = 0
            common_vector = []
            vector = vectors[mode].split()
            for res in nma_index:
                csum = 0
                for i in range(3):
                    ele = float(vector[res * 3 + i])
                    csum = csum+ele*ele
                    common_vector.append(ele)
                csum = sqrt(csum)
                correlationMode.append(csum)


            # Calculate the magnitude
            mag_mode = 0
            for r in common_vector:
                mag_mode += r * r

            mag_mode = sqrt(mag_mode)
            # Calculate Dot Product
            if len(common_vector) == len(delta_r):
                # print "Vectors Match"
                C = abs(np.corrcoef(correlationMode,correlationDR)[0,1])
                for i in range(len(common_vector)):
                    overlap += common_vector[i] * delta_r[i]

                overlap = overlap / (mag_d_r * mag_mode)
  
                spaces = len("mode: "+str(mode+1))
    
                spaces = 15-spaces
     
                if abs(overlap) in output:
                    output[abs(overlap)].append("Mode: " + str(mode+1) + ' '*spaces+ str(overlap) +'      '+str(C)+'\n')
                else:
                    output[abs(overlap)]=["Mode: " + str(mode+1) + ' '*spaces + str(overlap)+'      '+ str(C)+'\n']
            
      

    
        overlap_list = sorted(output.keys())
        overlap_list.reverse()


        w = open(outdir + "/" + outfile, 'w')
        w.write('MODE           Overlap              Correlation\n\n')
        for out in overlap_list:
            for o in output[out]:
                w.write(o)

        w.close()
    except IndexError:
        print ('\n**************************************\nFILE '+vtMatrix+' IS NOT A VALID EIGENVECTOR FILE:\n**************************************\n')
        sys.exit()



def calcDR(full_cords,empty_cords,count_common):

    struc = np.zeros((count_common,3))
    for i,res in enumerate(empty_cords):
        for j,c in enumerate(res):
            struc[i,j] = c

   
    target = np.zeros((count_common,3))
    for i,res in enumerate(full_cords):
        for j,c in enumerate(res):
            target[i,j] = c
    confAligned = sdrms.superpose3D(struc, target)
    alignedConf = confAligned[0]
    rmsd = confAligned[1]
    
    print ('\nRMSD between the two conformations = ' + str(rmsd) + '\n')

    # Calculate deltaR
    delta_r = []
    for j in range(len(alignedConf)):
        empty = alignedConf[j]
        full = target[j]
        rx = empty[0] - full[0]
        ry = empty[1] - full[1]
        rz = empty[2] - full[2]
        delta_r.append(rx)
        delta_r.append(ry)
        delta_r.append(rz)
    return delta_r

def commonCords(nma,lines_empty,common_residues):
    nma_index = []
    count_nma = 0
    count_common = 0
    empty_cords = []
    full_cords = []
    for line in lines_empty:
     
        info = line.split()
        atype = info[2].strip()
        res_type = info[3].strip()
        chain = info[4].strip()
        res = int(info[5].strip())
        if chain in common_residues:
            if res in common_residues[chain]:
                
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                cod = [x, y, z]
                empty_cords.append(cod)
                count_common+=1

    for line in nma:
   
        info = line.split()
        atype = info[2].strip()
        res_type = info[3].strip()
        chain = info[4].strip()
        res = int(info[5].strip())
        if chain in common_residues:
            if res in common_residues[chain]:
                    
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                cod = [x, y, z]
                full_cords.append(cod)
                nma_index.append(count_nma)
        count_nma+=1
    return(full_cords,nma_index,empty_cords,count_common)


def getCommonResidues(full_residues,empty_residues):
    common_residues = {}
    for ch in full_residues:
        for r in full_residues[ch]:
            if ch in empty_residues:
                if r in empty_residues[ch]:
                    if ch in common_residues:
                        common_residues[ch].append(r)
                    else:
                        common_residues[ch] = [r]
            else:
                break
    return common_residues


def parsePDB(pdb_file, atomT):

    #Open PDB FILE

    try:
        f = open(pdb_file, 'r')
        pdb_lines = f.readlines()
        f.close()
    except IOError:
        print ('\n**************************************\nFILE ' + pdb_file + ' NOT FOUND:\n**************************************\n')
        sys.exit()
  
    # determine the number of assymetric units
    number_of_protomers = 0
    currentResidue = 0
    currentChain = ''

    #sort residues by chain
    residues = {}

    #get lines of selected C atoms
    atom_lines = []

    for line in pdb_lines:
        if line.startswith("ATOM"):
            info = line.split()
            atype = info[2].strip()
            res_type = info[3].strip()
            chain = (info[4].strip())
            res = int(info[5].strip())
            if atype == atomT or (atype == "CA" and res_type == "GLY"):
                #Check Units
                if currentResidue == 0 and currentChain=='':
                    currentChain = chain
                    currentResidue = res
                    number_of_protomers+=1
                elif res==currentResidue and chain==currentChain:
                    number_of_protomers+=1

                #Check chains
                if chain in residues:
                    if res not in residues[chain]:
                        residues[chain].append(res)
                else:
                    residues[chain] = [res]
                #Append Atom Line
                atom_lines.append(line)
    return (number_of_protomers,residues,atom_lines)

    
def main(args):

    outdir = args.outdir
    outfile = args.output



    atomT = args.atomType.upper()
    if atomT!='CA' and atomT!='CB':
        print ('\n**************************************\nUnrecognised atom type\nInput Options:\nCA: to select alpha carbon atoms\nCB: to select beta carbon atoms\n**************************************')
        sys.exit()

    pdb_1 = args.pdbANM
    pdb_Conf = args.pdbConf
    vtMatrix = args.vtMatrix

    if pdb_1==pdb_Conf:
        print ('\n**************************************\nWARNING!!!\nConformational change PDB files are the same:\n--pdbANM: ' + pdb_1 + '\n--pdbConf: ' + pdb_Conf+ '\n**************************************\n')

    # determine the number of assymetric units and check compatability
    pdb1Info = parsePDB(pdb_1,atomT)
    pdbCInfo = parsePDB(pdb_Conf,atomT)
    
 
    number_of_protomersN = pdb1Info[0]
    full_residues = pdb1Info[1]
    nma = pdb1Info[2]
  
    number_of_protomersE = pdbCInfo[0]
    empty_residues =  pdbCInfo[1]
    lines_empty =  pdbCInfo[2]  

   
    if number_of_protomersN != number_of_protomersE:
        print ('\n**************************************\nERROR!!!:\nConformations contain a different number of assymetric units.\nCheck PDB files!!!\n**************************************\n')
        system.exit()


    common_residues = getCommonResidues(full_residues, empty_residues)

    commonK = sorted(common_residues.keys())

    totalC = 0
    for k in commonK:
        totalC +=len(common_residues[k])

    if not any(common_residues):
        print ('\n**************************************\nPDB Conformations are not compatable:\nSuggested error: Chain IDs do not match between PDB Files\n**************************************\n')
        system.exit()

   #Print warnings to user
  
    if len(full_residues.keys())>len(common_residues.keys()) or len(empty_residues.keys())>len(common_residues.keys()):
        print ('\n*****************************************************************\nWARNING!!!:\nNot all chains from PDB files were selected\nSuggested: Chain IDs do not match between PDB Files\n')
        print ("**************************************************************\nCorrelations calculated across "+str(totalC*number_of_protomersN)+" common residues ("+str(totalC)+" per "+str(number_of_protomersN)+" assymetric units).\nBreakdown per chain:\n")
        for k in commonK:
            print (k+": "+str(len(common_residues[k]))+ " residues per assymetric unit")
            resS=''
            for r in common_residues[k]:
                resS+=str(r)+" "
            print ("Residues selected include: "+resS+'\n')
        print("*****************************************************************\n")
    else:   
        print ("\n**************************************************************\nCorrelations calculated across "+str(totalC*number_of_protomersN)+" common residues ("+str(totalC)+" per "+str(number_of_protomersN)+" assymetric units).\nBreakdown per chain:\n")
        for k in commonK:
            print (k+": "+str(len(common_residues[k]))+ " residues per assymetric unit")
            resS=''
            for r in common_residues[k]:
                resS+=str(r)+" "
            print ("Residues selected include: "+resS+'\n')
        print("*****************************************************************\n")




    Common = commonCords(nma,lines_empty,common_residues)
    full_cords = Common[0]
    
    nma_index = Common[1]
    empty_cords = Common[2]
    count_common = Common[3]
    

    delta_r = calcDR(full_cords,empty_cords,count_common)

    calcConformation(delta_r,vtMatrix,nma_index,outdir,outfile)


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
    parser.add_argument("--pdbConf", help="")
    parser.add_argument("--pdbANM", help="")
    parser.add_argument("--vtMatrix", help="")  # note: change this from vtProtomer
    parser.add_argument("--output", help="Output file", default="ModesOfConformtionalChange.txt")
    parser.add_argument("--atomType", help="Enter CA to select alpha carbons or CB to select beta carbons", default='X')

    args = parser.parse_args()

    if args.welcome == "true":
        welcome_msg("Conformation mode", "Caroline Ross (caroross299@gmail.com)")

    print ('!=====================================================================================!')
    print ('! Please check the following:                                                         !')
    print ('! --pdbANM must be the PDB file that NMA was performed on                             !')
    print ('! --pdbConf must be a PDB of your complex aligned to the conformational change !')
    print ('!=====================================================================================!')

    # Check if required directories exist
    if not os.path.isdir(args.outdir):
        os.makedirs(args.outdir)

    # Check if args supplied by user
    if len(sys.argv) > 1:
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
