#!/usr/bin/env python

# visualiseVector.py
# Projects a set of eigenvectors of a normal onto the PDB structure. 
# Produces as set of frames to visualise the motion and Tcl script to plot the eigenvectors as a set of arrows
# Author: Caroline Ross: caroross299@gmail.com
# August 2017

import argparse
from datetime import datetime

from lib.utils import *


def main(args):

    atomT = args.atomType.upper()
    if atomT!='CA' and atomT!='CB':
	print '\n**************************************\nUnrecognised atom type\nInput Options:\nCA: to select alpha carbon atoms\nCB: to select beta carbon atoms\n**************************************'
	sys.exit()

    try:
        pdb_file = open(args.pdb, 'r')
        pdb_lines = pdb_file.readlines()
        pdb_file.close()
    except IOError:
        print '\n**************************************\nFILE '+args.pdb+' NOT FOUND:\n**************************************\n'
	sys.exit()

    # CHANGE THE FOLLOWING
    mode = args.mode

    structure = 'VISUAL'
 
    # get index of first atom
    for i in range(len(pdb_lines)):
        if pdb_lines[i].startswith("ATOM"):
            index_atom = i
            break
    header = pdb_lines[0:index_atom]
    all_atoms = pdb_lines[index_atom:]

    c_beta_atoms = []
    for atom in all_atoms:
        if atom.startswith("ATOM"):
            info = atom.split()
            first = info[0].strip()
            res = info[3]
            atom_type = info[2]
            if first == "ATOM" and (atom_type == atomT or (atom_type == "CA" and res == "GLY")):
                c_beta_atoms.append(atom)
        else:
            if "TER" in atom or "END" in atom:
                c_beta_atoms.append(atom)

    # Renumber the atoms
    for i in range(len(c_beta_atoms) - 1):
        a = c_beta_atoms[i]
        spaces = " " * (len(a[6:11]) - len(str(i + 1)))
        a = a[0:6] + spaces + str(i + 1) + a[11:]
        c_beta_atoms[i] = a

    # Determine Connections
    
    colours = ['red','blue','ochre','purple','yellow','red','cyan','pink','silver','violet','ochre','blue2','cyan2','iceblue','lime','green2','green3','violet','violet2','mauve']
    colourByChain = {}
    conect = []
    conect_chain = []
    atom1 = c_beta_atoms[0]
    chain1 = atom1.split()[4].strip()

    countColour = 0
    for i, atom in enumerate(c_beta_atoms, 1):
        if "TER" in atom or "END" in atom:
            continue
        chain = atom.split()[4].strip()
	if chain not in colourByChain and countColour<len(colours):
		colourByChain[chain]=colours[countColour]
		countColour+=1
        if chain == chain1:
            conect_chain.append(i)
        else:
            conect.append(conect_chain)
            conect_chain = []
            conect_chain.append(i)
            chain1 = chain
    conect.append(conect_chain)

    #Default to black if too many chains
    startColourBlack = False
    if countColour>len(colours):
	for chain in colourByChain:
            colourByChain[chain] = 'black'
	    startColourBlack = True

    keysC = colourByChain.keys()
    keysC.sort()
    print "ARROW COLOUR KEY BY CHAIN"
    for k in keysC:
        print k+': '+colourByChain[k]

    # Get the vectors
    # CHANGE HERE # may not be correct change, double check

    try:  
        vectorf = open(args.vectorFile, 'r')
        vectors = vectorf.readlines()
        vectorf.close()
    except IOError:
        print '\n**************************************\nFILE '+args.vectorFile+' NOT FOUND:\n**************************************\n'
	sys.exit()

    #Write VISUALISE
    try:
        w = open(args.outdir + "/" + "VISUALISE/" + structure + '_' + str(mode) + ".pdb", 'w')
        for i in range(0, 50):
            v_index = -1
            for atom in c_beta_atoms:
                if "ATOM" in atom:
                    v_index += 1 
                        
                    
                    v = vectors[v_index].split()
                    vx = float(v[0].strip())
                    vy = float(v[1].strip())
                    vz = float(v[2].strip())
                    x = round(float(atom[30:38].strip()) + (vx * i / 5), 3)
                    xspace = ' ' * (len(atom[30:38]) - len(str(x)))
                    y = round(float(atom[38:46].strip()) + (vy * i / 5), 3)
                    yspace = ' ' * (len(atom[38:46]) - len(str(y)))
                    z = round(float(atom[46:54].strip()) + (vz * i / 5), 3)
                    zspace = ' ' * (len(atom[46:54]) - len(str(z)))
                    atom = atom[0:30] + xspace + str(x) + yspace + str(y) + zspace + str(z) + atom[54:]
                    w.write(atom)
                else:
                    if "TER " in atom:
                        w.write(atom)
                    else:
                        if "END" in atom:
                            for Con in conect:
                                for c in range(len(Con) - 1):
                                    atom1 = str(Con[c])
                                    atom2 = str(Con[c + 1])
                                    w.write("CONECT" + " " * (5 - len(atom1)) +atom1 + " " * (5 - len(atom2)) + atom2 + "\n")
                            w.write(atom+'\n')

        w.close()

        # write arrows

        arrows = []
        steps = [0, 8]
        v_index = -1
        chainbreaks = []
        for cb in conect:
    	    chainbreaks.append(cb[-1])

        if startColourBlack:
	    arrows.append('proc vmd_draw_arrow {mol start end} {\n    set middle [vecadd $start [vecscale 0.9 [vecsub $end '
              '$start]]]\n    graphics $mol cylinder $start $middle radius 0.80\n    graphics $mol cone $middle $end '
              'radius 2.20\n}\ndraw color black\n')
        else:
            arrows.append('proc vmd_draw_arrow {mol start end} {\n    set middle [vecadd $start [vecscale 0.9 [vecsub $end '
              '$start]]]\n    graphics $mol cylinder $start $middle radius 0.80\n    graphics $mol cone $middle $end '
              'radius 2.20\n}\ndraw color '+colours[0]+'\n')

    

     
	nextAtom=[]
        for n,atom in enumerate(c_beta_atoms,1):
	    
            if "ATOM" in atom:
		chain = atom.split()[4].strip()
                v_index += 1
	        v = vectors[v_index].split()	
                vx = float(v[0].strip())
                vy = float(v[1].strip())
                vz = float(v[2].strip())
                x1 = round(float(atom[30:38].strip()) + (vx * steps[0]), 3)
                y1 = round(float(atom[38:46].strip()) + (vy * steps[0]), 3)
                z1 = round(float(atom[46:54].strip()) + (vz * steps[0]), 3)

                x2 = round(float(atom[30:38].strip()) + (vx * steps[1]), 3)
                y2 = round(float(atom[38:46].strip()) + (vy * steps[1]), 3)
                z2 = round(float(atom[46:54].strip()) + (vz * steps[1]), 3)
	
                arrows.append('draw arrow {' + str(x1) + ' ' + str(y1) + ' ' + str(z1) + '} {' + str(x2) + ' ' + str(y2) + ' ' + str(z2) + '}\n')
		if n in chainbreaks:
		    nextAtom = c_beta_atoms[n:]
		    for a in nextAtom:			
			if "ATOM" in a:
			    chain = a.split()[4].strip()
			    break		    
                    arrows.append('draw color '+colourByChain[chain]+'\n')
		    
  

        w = open(args.outdir + "/" + "VISUALISE/" + structure + '_ARROWS_' + str(mode) + ".txt", 'w')
        w.writelines(arrows)
        w.close()
    except IndexError:
        print '\n**************************************\nERROR!!\nVECTOR FILE and PDB FILE ARE NOT COMPATIBLE\n**************************************\n'
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
    # '3VBSFull_Aligned.pdb'
    parser.add_argument("--pdb", help="Coarse grained PDB file")  # '3VBSProtomer3_SCA.pdb'
    parser.add_argument("--mode", help="[int]", type=int)
    parser.add_argument("--vectorFile", help="File containing eigen vectors")  # 'ProtomerMode'
    parser.add_argument("--atomType", help="Enter CA to select alpha carbons or CB to select beta carbons", default='X')

    args = parser.parse_args()

    if args.welcome == "true":
        welcome_msg("Visualise vector", "Caroline Ross (caroross299@gmail.com)")

    # Check if required directories exist
    if not os.path.isdir(args.outdir):
        os.makedirs(args.outdir)
    if not os.path.isdir(args.outdir + "/" + "VISUALISE/"):
        os.makedirs(args.outdir + "/" + "VISUALISE/")

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
