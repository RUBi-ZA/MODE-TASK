#!/usr/bin/env python

# combinationMode.py
# Calculates the combined overlap and correlation for specified set of modes to a known conformational change
# Author: Caroline Ross: caroross299@gmail.com
# August 2017

import os
import sys
import argparse
from datetime import datetime
from utils import *
import numpy as np
from math import sqrt
import sdrms
import numpy as np

def main(args):

    atomT = args.atomType.upper()
    if atomT!='CA' and atomT!='CB':
	print '\n**************************************\nUnrecognised atom type\nInput Options:\nCA: to select alpha carbon atoms\nCB: to select beta carbon atoms\n**************************************'
	sys.exit()



    if args.pdbConf==args.pdbANM:
	print '\n**************************************\nWARNING!!!\nConformational change PDB files are the same:\n--pdbANM: '+args.pdbANM+'\n--pdbConf: '+args.pdbConf+'\n**************************************\n'
	

    try:
        f = open(args.pdbConf, 'r')
        lines_empty = f.readlines()
        f.close()
    except IOError:
        print '\n**************************************\nFILE '+args.pdbConf+' NOT FOUND:\n**************************************\n'
	sys.exit()



    try:
        f = open(args.pdbANM, 'r')
        nma = f.readlines()
        f.close()
    except IOError:
        print '\n**************************************\nFILE '+args.pdbANM+' NOT FOUND:\n**************************************\n'
	sys.exit()

    empty_residues = {}
    full_residues = {}
    common_residues = {}

    #determine the number of assymetric units
    number_of_protomersN = 0
    currentResidue = 0
    currentChain = ''
  

    for line in nma:
        if line.startswith("ATOM"):
            info = line.split()
            atype = info[2].strip()
            res_type = info[3].strip()
            chain = (info[4].strip())
            res = int(info[5].strip())
            if atype == atomT or (atype == "CA" and res_type == "GLY"):
		if currentResidue == 0 and currentChain=='':
			currentChain = chain
			currentResidue = res
			number_of_protomersN+=1	
		elif res==currentResidue and chain==currentChain:
			number_of_protomersN+=1
		
   
   

    #determine the number of assymetric units
    number_of_protomersE = 0
    currentResidue = 0
    currentChain = ''
  

    for line in lines_empty:
        if line.startswith("ATOM"):
            info = line.split()
            atype = info[2].strip()
            res_type = info[3].strip()
            chain = (info[4].strip())
            res = int(info[5].strip())
            if atype == atomT or (atype == "CA" and res_type == "GLY"):
		if currentResidue == 0 and currentChain=='':
			currentChain = chain
			currentResidue = res
			number_of_protomersE+=1	
		elif res==currentResidue and chain==currentChain:
			number_of_protomersE+=1
		
   
    if number_of_protomersN != number_of_protomersE:
        print '\n**************************************\nERROR!!!:\nConformations contain a different number of assymetric units.\nCheck PDB files!!!\n**************************************\n'
	system.exit()


   
    #determine co-ords of residues common in both PDB files		
    		
    for line in lines_empty:
        if line.startswith("ATOM"):
            info = line.split()
            atype = info[2].strip()
            res_type = info[3].strip()
            chain = info[4].strip()
            res = int(info[5].strip())
            if atype == atomT or (atype == "CA" and res_type == "GLY"):
		if chain in empty_residues:
                    if res not in empty_residues[chain]:
                        empty_residues[chain].append(res)
		    	
		else:
			empty_residues[chain] = [res]  
	
   	
    for line in nma:
        if line.startswith("ATOM"):
            info = line.split()
            atype = info[2].strip()
            res_type = info[3].strip()
            chain = info[4].strip()
            res = int(info[5].strip())
            if atype == atomT or (atype == "CA" and res_type == "GLY"):
	        if chain in full_residues:
                    if res not in full_residues[chain]:
                        full_residues[chain].append(res)
		   
		else:
		    full_residues[chain]=[res]

   

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

    commonK = common_residues.keys()
    commonK.sort()
    totalC = 0
    for k in commonK:
	totalC +=len(common_residues[k])




    if not any(common_residues):
        print '\n**************************************\nPDB Conformations are not compatable:\nSuggested error: Chain IDs do not match between PDB Files\n**************************************\n'
        system.exit()

   #Print warnings to user
  
    if len(full_residues.keys())>len(common_residues.keys()) or len(empty_residues.keys())>len(common_residues.keys()):
        print '\n*****************************************************************\nWARNING!!!:\nNot all chains from PDB files were selected\nSuggested: Chain IDs do not match between PDB Files\n'
        print "**************************************************************\nCorrelations calculated across "+str(totalC*number_of_protomersN)+" common residues ("+str(totalC)+" per "+str(number_of_protomersN)+" assymetric units).\nBreakdown per chain:\n"
        for k in commonK:
            print k+": "+str(len(common_residues[k]))+ " residues per assymetric unit"
	    resS=''
	    for r in common_residues[k]:
		    resS+=str(r)+" "
	    print "Residues selected include: "+resS+'\n'
        print"*****************************************************************\n"
    else:   
        print "\n**************************************************************\nCorrelations calculated across "+str(totalC*number_of_protomersN)+" common residues ("+str(totalC)+" per "+str(number_of_protomersN)+" assymetric units).\nBreakdown per chain:\n"
        for k in commonK:
            print k+": "+str(len(common_residues[k]))+ " residues per assymetric unit"
	    resS=''
	    for r in common_residues[k]:
		    resS+=str(r)+" "
	    print "Residues selected include: "+resS+'\n'
        print"*****************************************************************\n"


    count_common = 0
    empty_cords = []
    full_cords = []
    for line in lines_empty:
        if line.startswith("ATOM"):
            info = line.split()
            atype = info[2].strip()
            res_type = info[3].strip()
            chain = info[4].strip()
            res = int(info[5].strip())
	    if chain in common_residues:
		
                if res in common_residues[chain]:
                    if atype == atomT or (atype == "CA" and res_type == "GLY"):
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                        cod = [x, y, z]
                        empty_cords.append(cod)
			count_common+=1

    for line in nma:
        if line.startswith("ATOM"):
            info = line.split()
            atype = info[2].strip()
            res_type = info[3].strip()
            chain = info[4].strip()
            res = int(info[5].strip())
	    if chain in common_residues:		
                if res in common_residues[chain]:
                    if atype == atomT or (atype == "CA" and res_type == "GLY"):
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                        cod = [x, y, z]
                        full_cords.append(cod)

   
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
    
    print '\nRMSD between '+args.pdbANM+' and '+args.pdbConf+' = '+str(rmsd)+'\n'
   
    '''print alignedConf
    output = []
    print len(alignedConf)
    iatom = 0
    for atom in lines_empty:
	if atom.startswith('ATOM'):
            info = atom.split()
            atype = info[2].strip()
            res_type = info[3].strip()
            chain = info[4].strip()
            res = int(info[5].strip())
            if atype == atomT or (atype == "CA" and res_type == "GLY"):
                if chain in common_residues:
                    if res in common_residues[chain]:	
		        x = round(alignedConf[iatom][0],3)
                        xspace = ' ' * (len(atom[30:38]) - len(str(x)))
                        y = round(alignedConf[iatom][1],3)
                        yspace = ' ' * (len(atom[38:46]) - len(str(y)))
                        z = round(alignedConf[iatom][2],3)
                        zspace = ' ' * (len(atom[46:54]) - len(str(z)))
                        atom = atom[0:30] + xspace + str(x) + yspace + str(y) + zspace + str(z) + atom[54:]
                        output.append(atom)
		        iatom+=1
	else:
		output.append(atom)

    w = open("Test.pdb",'w')
    for out in output:
	w.write(out)
    w.close()'''
   

    
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

    # Get Eigenvectors for a mode
    # Calculate Indexes of vectors to be selected
    interface_index = []
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
                        interface_index.append(count)
                count += 1


    try:   
        f = open(args.vtMatrix, 'r')
        vectors = f.readlines()
        f.close()
    except IOError:
        print '\n**************************************\nFILE '+args.vtMatrix+' NOT FOUND:\n**************************************\n'
	sys.exit()

    directions = []
    mode_range = []
    modeList = args.modes
    if modeList =='0':
	print '\n**************************************\nERROR!!!\nNo Modes Specified\nEnter mode numbers separated by a comma\nEg: --modes 1,5,7\n**************************************\n'
    else:
	try:
	    modesL = modeList.split(',')
	    for m in modesL:
		mode_range.append(int(m)-1)
	    
	except ValueError:
	    print '\n**************************************\nERROR!!!\nUnrecognised input\nEnter mode numbers separated by a comma\nEg: --modes 1,5,7\n**************************************\n'
	    
    
	
    
    output = {}
    try:
        for mode in mode_range:
	    correlationMode = []
            overlap = 0
            common_vector = []
            vector = vectors[mode].split()
            for res in interface_index:
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
		
	    	if overlap <0:
		    directions.append(-1)
		else:
		    directions.append(1)
                spaces = len("mode: "+str(mode+1))
	    
	        spaces = 15-spaces
	    
	        if abs(overlap) in output:
                    output[abs(overlap)].append("Mode: " + str(mode+1) + ' '*spaces+ str(overlap) +'      '+str(C)+'\n')
	        else:
		    output[abs(overlap)]=["Mode: " + str(mode+1) + ' '*spaces + str(overlap)+'      '+ str(C)+'\n']

	    else:
		print '\n**************************************\nConformation PDB and ANM PDB are not compatable\nExiting...\n**************************************\n'
            
    	    

    
        overlap_list = output.keys()
        overlap_list.sort()
        overlap_list.reverse()



	#Calculate the resultant vector
	OverlapC = 0
        correlationC = 0
	correlationMode = []
	vectorR = []
	for res in interface_index:
	    vectorR.append(0)
	    vectorR.append(0)
	    vectorR.append(0)

    	for d,mode in enumerate(mode_range):
	    direct=directions[d]
	    vector = vectors[mode].split()
	    CommonVector = []
	    for res in interface_index:
               
		for i in range(3):
		    ele = direct*float(vector[res*3+i])
 		  
		    CommonVector.append(ele)
	    
	        
	    for j,v in enumerate(CommonVector):
		vectorR[j]+=v

        mag_mode = 0
        for r in vectorR:
                mag_mode += r * r
	mag_mode = sqrt(mag_mode)

	
	for c in range(0,len(vectorR),3):
	    csum = vectorR[c]*vectorR[c] + vectorR[c+1]*vectorR[c+1] +vectorR[c+2]*vectorR[c+2]
	    csum = sqrt(csum)
	    correlationMode.append(csum)

	print "CHECKING "+str(len(correlationMode))+"\n*******\n"
	if len(vectorR)==len(delta_r):
	  
	    for i in range(len(vectorR)):
		OverlapC+=vectorR[i]*delta_r[i]
		
	    OverlapC = OverlapC/(mag_d_r * mag_mode)
	    correlationC = abs(np.corrcoef(correlationMode,correlationDR)[0,1])

	    print OverlapC
	    print  correlationC
	else:
	    print '\n**************************************\nConformation PDB and ANM PDB are not compatable\nExiting...\n**************************************\n'



   



        w = open(args.outdir + "/" + args.output, 'w')
        w.write('MODE           Overlap              Correlation\n\n')
        for out in overlap_list:
	    for o in output[out]:
                w.write(o)
	w.write('\n********************************************************************************************************\nCombinded Overlap = '+str(OverlapC)+'\nCombinded Correlation = '+str(correlationC)+'\n********************************************************************************************************\n')

        w.close()
    except IndexError:
        print '\n**************************************\nFILE '+args.vtMatrix+' IS NOT A VALID EIGENVECTOR FILE:\n**************************************\n'
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
    parser.add_argument("--log-file", help="Output log file (default: standard output)", default=None)
    parser.add_argument("--outdir", help="Output directory", default="output")

    # custom arguments
    parser.add_argument("--pdbConf", help="")
    parser.add_argument("--pdbANM", help="")
    parser.add_argument("--vtMatrix", help="")  # note: change this from vtProtomer
    parser.add_argument("--output", help="Output file", default="ModeSpecificConformationalChange.txt")
    parser.add_argument("--atomType", help="Enter CA to select alpha carbons or CB to select beta carbons", default='X')
    parser.add_argument("--modes", help="Calculate the overlap for a combination of specific modes\nEnter mode numbers separated by a comma\nEg: --modes 1,5,7", default='0')

    args = parser.parse_args()

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
        print "No arguments provided. Use -h to view help"
