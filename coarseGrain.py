#!/usr/bin/env python

# coarseGrain.py
# Selects a sub-set of carbon atoms, equally distributed across the complex
# Author: Caroline Ross: caroross299@gmail.com
# August 2017


import os
import sys
import argparse
from datetime import datetime

from lib.utils import *

import math


def parseAssembly(pdb_file,atomT):
    try:
        f = open(pdb_file, 'r')
        lines = f.readlines()
        f.close()
    except IOError:
        print ('\n**************************************\nFILE '+pdb_file+' NOT FOUND\n**************************************\n')
        sys.exit()

# get index of first atom
    for i in range(len(lines)):
        if lines[i].startswith("ATOM"):
            index_atom = i
            break
    header = lines[0:index_atom]
    all_atoms = lines[index_atom:]

    chain_dics = {}
    number_protomer_atoms = 0
    c_atoms = []
    c_cords = []

    for atom in all_atoms:
        if atom.startswith("ATOM"):
            coords = [] 
            info = atom.split()
            chain = info[4]
            res = info[3]
            atom_type = info[2]
            if atom_type == atomT or (atom_type == "CA" and res == "GLY"): #Give user the option for C-aplhas or C-betas
                if number_protomer_atoms == 0:
                    chain1 = chain
                if chain not in chain_dics:
                    number_protomer_atoms += 1
                if chain != chain1:
                    if chain1 not in chain_dics:
                        chain_dics[chain1] = 1
                    else:
                        # Is macro molecule
                        chain_dics[chain1] += 1
                    chain1 = chain
                c_atoms.append(atom)
                coords.append(float(atom[30:38].strip().strip()))
                coords.append(float(atom[38:46].strip().strip()))
                coords.append(float(atom[46:54].strip().strip()))
                c_cords.append(coords)
        if atom.startswith("END") and not "ENDMDL" in atom:
            if chain1 in chain_dics:
                chain_dics[chain1] += 1
            else:
                chain_dics[chain1] = 1
    number_of_protomers = chain_dics[chain1]

    return (header,number_protomer_atoms,number_of_protomers,c_atoms,c_cords)

def coarseGrain(c_g,starting_atom,number_protomer_atoms, number_of_protomers, c_cords,c_atoms):


    index_of_selected_atoms = []
    starting_atom_i = starting_atom - 1  # Index for starting atom
    
    
    protomer_c_betas = c_cords[0:number_protomer_atoms]
    coords_start = protomer_c_betas[starting_atom_i]
    distances_from_start = []
    distance_index = {}  # holds the index of the atoms in order of distance from

    xstart = coords_start[0]
    ystart = coords_start[1]
    zstart = coords_start[2]

    for i in range(len(protomer_c_betas)):
        atom = protomer_c_betas[i]
        if i == starting_atom_i:
            continue
        x = atom[0]
        y = atom[1]
        z = atom[2]
        distance = ((xstart - x) * (xstart - x)) + ((ystart - y)
                                                    * (ystart - y)) + ((zstart - z) * (zstart - z))
        distance = math.sqrt(distance)
        distances_from_start.append(distance)
        distance_index[distance] = i
    distances_from_start.sort()


   
    index_of_selected_atoms = []
    index_of_selected_atoms.append(starting_atom_i)

    if c_g == 1:
        c_g_select = 0
    elif c_g ==2:
        c_g_select = 1
    else:
        c_g_select = (c_g * (c_g - 1)) - c_g
    
    # selects atoms which are not within this distance to already selected atoms
    try:
        cutoff = distances_from_start[c_g_select]
        atom_index = distance_index[cutoff]
        index_of_selected_atoms.append(atom_index)
    
    
        distribution = {}
    
        #loops through ordered list and selects all suitable atoms
        for dist in distances_from_start[c_g_select+ 1:]:
            atom_index = distance_index[dist]
            x = protomer_c_betas[atom_index][0]
            y = protomer_c_betas[atom_index][1]
            z = protomer_c_betas[atom_index][2]
            too_close = False
            local_distribution = []
            for atom in index_of_selected_atoms:
                x1 = protomer_c_betas[atom][0]
                y1 = protomer_c_betas[atom][1]
                z1 = protomer_c_betas[atom][2]
                surrounding_distance = ((x1 - x) * (x1 - x)) + ((y1 - y) * (y1 - y)) + ((z1 - z) * (z1 - z))
                surrounding_distance = math.sqrt(surrounding_distance)
                local_distribution.append(surrounding_distance)
                if surrounding_distance < cutoff:
                    too_close = True
                    break
            if not too_close:
                index_of_selected_atoms.append(atom_index)
                local_distribution.sort()
                distribution[atom_index] = local_distribution
    
            # print index_of_selected_atoms
        print ("------------------------------------------------------------")
        print ("SUMMARY OF COARSE GRAINING PERFORMED AT LEVEL "+str(c_g))
        print ("No. atoms selected per unit: " + str(len(index_of_selected_atoms)) + " from " + str(number_protomer_atoms) + " orignal residues")
        print ("No. atoms selected per macromolecule: " + str(len(index_of_selected_atoms) * number_of_protomers) + " from " + str(number_protomer_atoms* number_of_protomers) + " orignal residues")
        print ("------------------------------------------------------------")
    

    
        index_of_selected_atoms.sort()
        # write a pdb file lines of the Coarse-Grained Capsid with renumbered atoms
        selected_c_beta_lines = []
        # Includes all c atoms of the first pentamer and then coarse grains the rest of the surrounding capsid
        count = 0
    
        for i in range(0, number_of_protomers):  # parameter
            for j in index_of_selected_atoms:
                index = j + i * number_protomer_atoms
                a = c_atoms[index]
                spaces = " " * (len(a[6:11]) - len(str(count + 1)))
                a = a[0:6] + spaces + str(count + 1) + a[11:]
                selected_c_beta_lines.append(a)
                count += 1
            ter = "TER" + " " * 3 + " " * (5 - len(str(count))) + str(count) + " " * 6 + a.split()[3] + " " + a.split()[4] + " " + " " * (3 - len(a.split()[5])) + a.split()[5] + " \n"
            selected_c_beta_lines.append(ter)
        return selected_c_beta_lines
    except IndexError:
       print ("\n**************************************\nERROR: Coarse Grain Level = "+str(c_g)+" is TOO LARGE\nNo output generated\n**************************************")
        

def writeCG(outfile,outdir,header,selected_c_beta_lines,c_g):
 
    if outfile == "ComplexCG.pdb":
        outfile = outfile.strip('.pdb')+str(c_g)+'.pdb'     
    elif '.pdb' in outfile:
        outfile = outfile.strip('.pdb')+'CG'+str(c_g)+'.pdb'    
	    
    elif not '.pdb' in outfile and not'.' in outfile:
        outfile+='CG'+str(c_g)
        outfile+='.pdb'
    elif not '.pdb' in outfile and '.' in outfile:
        print ('\n**************************************\nSpecified output file name is not comptable with PDB format\nDefault file created: ComplexCG_'+str(c_g)+'.pdb\n**************************************')
        outfile = 'ComplexCG'+str(c_g)+'.pdb'
	
    
    w = open(args.outdir + "/" + outfile, 'w')
    w.writelines(header)
    w.writelines(selected_c_beta_lines)
    w.write("END")
    w.close()


def main(args):


    #Get Input Parameters
    pdb_file = args.pdb
    atomT = args.atomType.upper()
    if atomT!='CA' and atomT!='CB':
        print ("\n**************************************\nUnrecognised atom type\nInput Options:\nCA: to select alpha carbon atoms\nCB: to select beta carbon atoms\n**************************************\n")
        sys.exit()


    # Set level of coarsegrain
    c_gList = []
    try:
        cgs = args.cg.split(',')
        for c_g in cgs:
            c_gList.append(int(c_g.strip()))
    except ValueError:
        print ("\n**************************************\nInput for Coarse Graining Level: "+str(args.cg)+" is INVALID\nDefault CG Level = 4 has been used\n**************************************\n")
        c_gList =[]
        c_gList.append(4)

    starting_atom = args.startingAtom  # residue number of starting atoms
    if starting_atom<=0:        
        print ("\n**************************************\nStarting Atom: "+str(starting_atom)+" is INVALID\nDefault Starting Atom = 1 has been used\n**************************************\n")
        starting_atom = 1
    
    outfile = args.output
    outdir = args.outdir
    #Parse PDB File      
    cAssembly = parseAssembly(pdb_file,atomT)
    header  = cAssembly[0]
    number_protomer_atoms = cAssembly[1]
    number_of_protomers = cAssembly[2]
    c_atoms = cAssembly[3]
    c_cords = cAssembly[4]
    #Perform Coarse Graining
    for c_g in c_gList:
        if c_g == 0 or c_g<0:
            print ("\n**************************************\nERROR: Coarse Grain Level = "+str(c_g)+" is INVALID\nNo output generated\n**************************************")
            continue
        selected_c_beta_lines = coarseGrain(c_g,starting_atom,number_protomer_atoms, number_of_protomers, c_cords,c_atoms) 
    
        if selected_c_beta_lines:
            writeCG(outfile,outdir,header,selected_c_beta_lines,c_g)
    

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

    # standard arguments for logging and output
    parser.add_argument("--silent", help="Turn off logging", action='store_true', default=False)
    parser.add_argument("--welcome", help="Display welcome message (true/false)", default="true")
    parser.add_argument("--log-file", help="Output log file (default: standard output)", default=None)
    parser.add_argument("--outdir", help="Output directory", default="output")

    # custom arguments
    parser.add_argument("--output", help="File name for Coarse Grained PDB", default='ComplexCG.pdb')
    parser.add_argument("--pdb", help="PDB input file")
    parser.add_argument("--cg", help="Course grain levels: Increase level to increase amount of coarse graining:\nEnter either a single level e.g. 4 or a comma separted list for the production of multiple CG models e.g 3,4,5", default='4')
    parser.add_argument("--startingAtom", help="Residue number of first carbon atom to be selected [int]", default=1, type=int)
    parser.add_argument("--atomType", help="Enter CA to select alpha carbons or CB to select beta carbons", default='CA')
    # parser.add_argument("--protomerAtoms", help="", default=0, type=int)

    args = parser.parse_args()

    if args.welcome == "true":
        welcome_msg("Coarse grain", "Caroline Ross (caroross299@gmail.com)")

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
