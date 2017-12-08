#!/usr/bin/env python

# assemblyCovariance.py
# Calcualtes covariance matrices for the following:
# 1) over all modes
# OR
# 2) specified modes
# AND
# 3) For a single assymteric unit
# OR
# 4) For a cluster of specified assymetric units
# Author: Caroline Ross: caroross299@gmail.com
# August 2017

# Calculates and Returns Diagonals of Correlated Matrix for a given set of modes

import os
import argparse
from datetime import datetime

from lib.utils import *
import numpy as np

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm as CM
import numpy as np


def plotCorrelation(M,fname,Vn,Vx,output):
    M = np.array(M) #converts M to and np array. matplotlib requires it to be in this format
    fig, ax = plt.subplots() #creates a new empty figure
    colors = [('white')] + [(CM.jet(i)) for i in xrange(40,250)] #defines the color scheme of the map
    new_map = matplotlib.colors.LinearSegmentedColormap.from_list('new_map', colors, N=300) #define the type of map to be plotted - a LinearSegmentedColormap
    heatmap = ax.pcolor(M, cmap=new_map,vmin=-Vn, vmax=Vx) #the next two lines create the map based on M and the colors, vmin and vmax set the range of values plotted
    heatmap = plt.pcolor(M, vmin=-Vn, vmax=Vx)
    fig = plt.gcf() #sets up figure
    fig.set_size_inches(10,15) #sets sixe of figure and the next few lines determines properties about the axis
    ax.set_frame_on(False)
    ax.invert_yaxis()
    ax.xaxis.tick_top()
    # rotate the
    #plt.xticks(rotation=90)
    ax.grid(False)
    # Turn off all the ticks
    ax = plt.gca()
    for t in ax.xaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False
    for t in ax.yaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False

    cbar = plt.colorbar(heatmap,orientation="horizontal") #actually plots the heatmap figure which has been created
    plt.savefig(output+"/"+fname,dpi=300) #save fig
    plt.close('all') #close fig    

def calcCovariance(wMatrix, vtMatrix, modes, definedUnits, specifiedUnits, zoom, vmin, vmax,output):  

    try: 
        f = open(wMatrix, 'r')
        eigen_values = f.readlines()
        f.close()
    except IOError:
        print ('\n**************************************\nERROR!! '+protein_name+' W-MATRIX FILE: '+wMatrix+' NOT FOUND:\n**************************************\n')
        sys.exit()
    try: 
        f = open(vtMatrix, 'r')
        eigen_vectors = f.readlines()
        f.close()
    except IOError:
        print ('\n**************************************\nERROR!! '+protein_name+' VT-MATRIX FILE: '+vtMatrix+' NOT FOUND:\n**************************************\n')
        sys.exit()

    total_modes = len(eigen_values)
    w_inv = np.zeros((total_modes, total_modes))

    try:

        if modes == 'all':
        
    
            #Create A Full W Inverse Matrix (This is if we want correlation averaged over all modes)
        
            for i in range(total_modes):
                if i >5:
                    w_inv[i, i] = 1 / (float(eigen_values[i].split()[1].strip()))
        else:
            for i in modes:
                w_inv[i, i] = 1 / (float(eigen_values[i].split()[1].strip()))
    except TypeError:
        print ('\n**************************************\nERROR!! INVALID INPUT FOR MODE RANGE: DEFAULT all Has been selected\n**************************************\n')
        w_inv = np.zeros((total_modes, total_modes))
        for i in range(total_modes):
            if i >5:
                w_inv[i, i] = 1 / (float(eigen_values[i].split()[1].strip()))

    # Read In U and VT full Matrix as U is the transpose of VT I only read in VT and create U from the VT matrix
    # info. So we can exclude U output from C++ script for faster analysis
    print ("Parsing VT_Matrix")

    try:
        fvt = open(vtMatrix, 'r')  
        eigen_vectors = fvt.readlines()
        fvt.close()
    except IOError:
        print ('\n**************************************\nERROR!! '+protein_name+' VT-MATRIX FILE: '+vtMatrix+' NOT FOUND:\n**************************************\n')
        sys.exit()

    print ("Calculating Transpose of VT")
    v_t = np.zeros((total_modes, total_modes))
    u = np.zeros((total_modes, total_modes))

    for i in range(total_modes):
        vectors = eigen_vectors[i].split()
        for j in range(total_modes):
            vector = float(vectors[j].strip())
            v_t[i, j] = vector
            u[j, i] = vector

    # Calculate Correlation Matrices
    print ("Calculating Covariance Matrix")
    w_v_t = np.dot(w_inv, v_t)
    c = np.dot(u, w_v_t)

    # Calculate Trace of the Correlation Matrices
    trace_c = np.zeros((total_modes / 3, total_modes / 3))
  

    for i in range(0, total_modes, 3):
        for j in range(0, total_modes, 3):
            trace = 0
            for k in range(3):
                trace = trace + c[i + k, j + k]
            trace_c[i / 3, j / 3] = trace

    print ("Covariance Matrix Calculation Complete")

    #Construct full image
    plotCorrelation(trace_c,"Covariance_Full.png",vmin,vmax,output)
    w=open(output+'/Covariance_Full.txt','w')
    for i in range(0, total_modes/3):
        for j in range(0, total_modes/3):
            w.write(str(trace_c[i,j])+" ")
        w.write('\n')
    w.close()


    try:
        #Construct assymteric image
        numberUnits = len(specifiedUnits)
        sizePerUnit = len(definedUnits[1])

        trace_cf = np.zeros((numberUnits*sizePerUnit,numberUnits*sizePerUnit))
        cf_rows = 0
        cf_cols = 0
    
        for i in specifiedUnits:
            startIndex = (i-1)*sizePerUnit
            endIndex = i*sizePerUnit
            for k in range(startIndex,endIndex):
                cf_cols = 0
                for j in specifiedUnits:
                    startIndexC = (j-1)*sizePerUnit
                    endIndexC = j*sizePerUnit
                    for l in range (startIndexC,endIndexC):
                        trace_cf[cf_rows,cf_cols] = trace_c[k,l]
                        cf_cols+=1
                cf_rows+=1 
        plotCorrelation(trace_cf,"Covariance_AUnits.png",vmin,vmax,output)
        w=open(output+'/Covariance_AUnits.txt','w')
        for i in range(numberUnits*sizePerUnit):
            for j in range(numberUnits*sizePerUnit):
                w.write(str(trace_cf[i,j])+" ")
            w.write('\n')
        w.close()
        
  
        if zoom!='0,0':
            unit = int(zoom.split(',')[0])
            chain = zoom.split(',')[1]
            unitAtoms = definedUnits[unit]

            count = 0
            start = 0
            end = 0
            for atom in unitAtoms:
                count+=1
                ch= atom.split()[4]

                if ch == chain and start==0:
                    start = count    
                elif ch!=chain and start!=0:
                    end = count
                    break
                elif atom == unitAtoms[-1]:
                    end = count

            start = start-1+ (unit-1)*sizePerUnit
            end = end + (unit-1)*sizePerUnit
            trace_Z = trace_c[start:end,start:end] 
            plotCorrelation(trace_Z,"Covariance_Zoom.png",vmin,vmax,output)
            w=open(output+'/Covariance_Zoom.txt','w')
            for i in range(end-start):
                for j in range(end-start):
                    w.write(str(trace_Z[i,j])+" ")
                w.write('\n')
            w.close()
    except IndexError:
        print ('\n**************************************\nERORR!! PDB FILE AND ANM MATRICES ARE IMCOMPATABLE\nCHECK INPUT PARAMETERS FOR:\n1) INCORRECT PDB FILE \n2) INCORRECT MATRICES \n3) INCORRECT SPECIFIED ASSYMETRIC UNITS\n**************************************\n')
        sys.exit()            



def defineAssemetricUnits(pdb,atom):

    try:
        pdb_file = open(pdb, 'r')
        pdb_lines = pdb_file.readlines()
        pdb_file.close()
    except IOError:
        print ('\n**************************************\nERROR!! FILE ' +pdb + ' NOT FOUND:\n**************************************\n')
        sys.exit()


    protomers = {}
    
    #Determine protomers
    protomerLines = []
    chain_dics = []
    number_protomer_atoms = 0
    countProtomers = 1
    for l in pdb_lines:

        if l.startswith("ATOM"):
            info = l.split()
            start = info[0]
            res_type = info[3].strip()
            aT = info[2]
            chain = l.split()[4]
            if aT == atom or (aT == "CA" and res_type == "GLY"):       
                if number_protomer_atoms == 0:
                    chain1 = chain
                if chain not in chain_dics:
                    protomerLines.append(l)
                    number_protomer_atoms+=1
                    #print protomerLines
                if chain != chain1:
                    chain_dics.append(chain1)
                    chain1 = chain
                if chain in chain_dics:
                    #Is macro molecule
                    protomers[countProtomers] = protomerLines
                    countProtomers+=1
                    chain_dics = []
                    protomerLines = []
                    protomerLines.append(l)
                    #chain1 = chain
  
            if l.startswith("TER"):
                protomerLines.append(l)
            if l.startswith('END'):
               protomerLines.append(l)
    
        protomers[countProtomers] = protomerLines
    return protomers



def main(args):

    pdbFile = args.pdb
    modes = args.modes
    WM = args.wMatrix
    VTM = args.vtMatrix
    AUnits = args.assymetricUnits
    zoomIndex = args.zoom
    aType = args.atomType
    vmin = args.vmin
    vmax = args.vmax
    output = args.outdir

    #Check paramaters
    atomT = args.atomType.upper()
    if atomT!='CA' and atomT!='CB':
        print ('\n**************************************\nERROR!! Unrecognised atom type\nInput Options:\nCA: to select alpha carbon atoms\nCB: to select beta carbon atoms\n**************************************')
        sys.exit()



    modes = args.modes
    mode_range = []
    try:
        if modes !='all':
    
            if ':' in modes:
                mRange = modes.split(':')
                mode_range = range(int(mRange[0])-1,int(mRange[1]))
            else:
                modesSpecific = modes.split(',')
                for m in modesSpecific:
                    if int(m) <= 0 :
                        print ('\n**************************************\nWARNING!! ZERO OR NEGATIVE MODES SPECIFIED. MATRIX INDEXES WILL BE AFFECTED\nRESULTS MAY NOT BE ACCURATE\n**************************************')
                    mode_range.append(int(m)-1)
        else:
            mode_range = 'all'
    except TypeError:
        print ('\n**************************************\nERROR!! INVALID INPUT FOR MODES: Default all has been selected\n**************************************')
        mode_range = 'all'

    AUnitsL = []
    try:
        if AUnits == '1':
            AUnitsL.append(1)
        else:
            AUnits = AUnits.split(',')
            for a in AUnits:
                if int(a) <= 0 :
                        print ('\n**************************************\nWARNING!! ZERO OR NEGATIVE ASSYMETRIC UNITS. MATRIX INDEXES WILL BE AFFECTED\nRESULTS MAY NOT BE ACCURATE\n**************************************')
                AUnitsL.append(int(a))
    except TypeError:
        print ('\n**************************************\nERROR!! INVALID INPUT FOR ASSYMETRIC UNITS: Default 1 has been selected\n**************************************')
        AUnitsL = []
        AUnitL.append(1)

    definedAUnits = defineAssemetricUnits(pdbFile,aType)
    calcCovariance(WM, VTM, mode_range, definedAUnits, AUnitsL, zoomIndex,vmin,vmax,output)






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
    parser.add_argument("--pdb", help="Input PDB file that was analysed by ANM") 
    parser.add_argument("--modes", help="1) Calculate the covariance matrix over all modes by using the option all:\n  E.g --modes all\nOR\n2) Enter a select range of modes in format M1:M2\n  E.g To calculate the covariance matrix over the first 20 non-zero modes enter --modes 7:27\nOR\n3) Calculate the covariance matrix for a combination of specific modes\nEnter mode numbers separated by a comma\n  E.g: --modes 1,5,7", default="all") 
    parser.add_argument("--wMatrix", help="Text file of Eigevalues of pdb, in format output from ANM.cpp")
    parser.add_argument("--vtMatrix", help="Text file of Eigevectors of pdb, in row (VT) format output from ANM.cpp")
    parser.add_argument("--assymetricUnits", help="1) Specify a single unit.\n  E.g --assymentricUnit 1 \nOR\n 2) Provide a list of assymteric units.\n  E.g --assymetricUnitsText 1,4,5", default = '1')
    parser.add_argument("--zoom", help="1) Zoom into a single chain within a specified assymetric unit.\nSpecify the unit and then the chain in a comma separted list single unit.\n  E.g --zoom 1,2 will zoom into chain 2 in assymteric unit 1",default = '0,0')
    parser.add_argument("--atomType", help="Enter CA to select alpha carbons or CB to select beta carbons", default='CA')
    parser.add_argument("--vmin", help="[int] minimum axes value", type=int, default = -0.1)
    parser.add_argument("--vmax", help="[int] maximum axes value", type=int, default = 0.1)
    args = parser.parse_args()

    if args.welcome == "true":
        welcome_msg("Assembly Covariance", "Caroline Ross (caroross299@gmail.com)")

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

