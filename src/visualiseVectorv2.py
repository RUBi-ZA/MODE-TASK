#!/usr/bin/env python
"""
Projects a set of eigenvectors of a normal onto the PDB structure. Offers the
user additional options for visualisation. Optionally produces a set of frames
to visualise the motion, Tcl and NGL scripts to plot the eigenvectors as a set
of arrows in VMD and NGLView respectively

Authors: Caroline Ross: caroross299@gmail.com
(with some modifications from Olivier Sheik Amamuddy)
18th Feb 2021
"""
import math
import argparse
from datetime import datetime
from lib.utils import *
from visualisevector2ngl import main as visualisevector2ngl

class Args:
    """Dummy class for handling filenames"""
    def __init__(self):
        self.filename = None
        self.outfilename = None

def get_nmodes(vtMatrix):
    """Returns the number of modes from the
    matrix of eigenvectors"""
    with open(vtMatrix, "r") as f:
        nmodes = f.readline().rstrip().split()
    return len(nmodes)

def getVector(vtMatrix, mode, direction):
    """Selects a row (eigenvector) from the VT matrix and
    reconstitutes each of the displacements vectors"""
    try:
        f = open(vtMatrix, "r")
        vlines = f.readlines()
        f.close()
    except IOError:
        print("ERROR: {} NOT FOUND:".format(vtMatrix))
        sys.exit(1)
    modeVec = vlines[mode].split()
    vectorFile = []
    numVec = len(modeVec)
    for i in list(range(0, numVec, 3)):
        vec = []
        for j in list(range(3)):
            vec.append(float(modeVec[i+j].strip()))
        vectorFile.append(vec)
    unitVectorFile = []
    for v in vectorFile:
        unitVec = []
        mag = 0
        for e in v:
            mag = mag+(e*e)
        mag = math.sqrt(mag)
        for e in v:
            e = direction*(e/mag)
            unitVec.append(e)
        unitVectorFile.append(unitVec)
    return unitVectorFile

def defineAsymmetricUnits(pdb_lines, vectors, aUnits, atom):
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
        protomers[countProtomers] = protomerLines

    Cvectors = []
    c_atoms = []
    atoms_per_protomer = len(protomers[1])

    for unit in aUnits:
        unit = int(unit)
        unit_atoms = protomers[unit]
        for i,a in enumerate(unit_atoms):
            c_atoms.append(a)
            Cvectors.append(vectors[i+(atoms_per_protomer*(unit-1))])
    c_atoms.append("END")
    return (c_atoms,Cvectors)

def visualise(pdb_lines, atomT, mode, vectors, outdir, head, tail,
              arrowLength, coloursSpecified, chainSpecified, aUnits,
              animate=False):
    args = Args()
    if aUnits[0] == 'all':
        structure = 'VISUAL'
        # get index of first atom
        for i in range(len(pdb_lines)):
            if pdb_lines[i].startswith("ATOM"):
                index_atom = i
                break
        #header = pdb_lines[0:index_atom]
        all_atoms = pdb_lines[index_atom:]
        c_beta_atoms = []
        for atom in all_atoms:
            if atom.startswith("ATOM"):
                info = atom.split()
                first = info[0].strip()
                res = info[3]
                atom_type = info[2]
                if first == "ATOM" and (atom_type == atomT or \
                        (atom_type == "CA" and res == "GLY")):
                    c_beta_atoms.append(atom)
            else:
                if "TER" in atom or "END" in atom:
                    c_beta_atoms.append(atom)
    else:
        structure = 'VISUAL_AUNITS'
        new_info = defineAsymmetricUnits(pdb_lines, vectors, aUnits, atomT)
        c_beta_atoms = new_info[0]
        vectors = new_info[1]

    # Renumber the atoms
    for i in range(len(c_beta_atoms) - 1):
        a = c_beta_atoms[i]
        spaces = " " * (len(a[6:11]) - len(str(i + 1)))
        a = a[0:6] + spaces + str(i + 1) + a[11:]
        c_beta_atoms[i] = a

    # Determine Connections
    if coloursSpecified == 'none':
        colours = ['red', 'blue', 'ochre', 'purple', 'yellow', 'red', 'cyan',
                   'pink', 'silver', 'violet', 'ochre', 'blue2', 'cyan2',
                   'iceblue', 'lime', 'green2', 'green3', 'violet', 'violet2',
                   'mauve']
    else:
        colours = coloursSpecified.strip().split(',')
        for i,col in enumerate(colours):
            col = col.lower()
            colours[i] = col
    colourByChain = {}

    if chainSpecified !='all':
        colourByChain[chainSpecified]=colours[0]

    conect = []
    conect_chain = []
    atom1 = c_beta_atoms[0]
    chain1 = atom1.split()[4].strip()
    startColourBlack = False
    countColour = 1
    for i, atom in enumerate(c_beta_atoms, 1):
        if "TER" in atom or "END" in atom:
            continue
        chain = atom.split()[4].strip()
        if chain not in colourByChain and \
                countColour<=len(colours) and \
                chainSpecified =='all':
            colourByChain[chain]=colours[countColour-1]
            countColour+=1
        elif chain not in colourByChain and \
                countColour>len(colours) and \
                chainSpecified =='all':
            startColourBlack = True
        if chain == chain1:
            conect_chain.append(i)
        else:
            conect.append(conect_chain)
            conect_chain = []
            conect_chain.append(i)
            chain1 = chain
    conect.append(conect_chain)
    #Default to black if too many chains
    if startColourBlack:
        print('WARNING: User may have specified less colours than number of chains.\
               All arrows are now coloured black')
        colourByChain = {}
        for i, atom in enumerate(c_beta_atoms, 1):
            if "TER" in atom or "END" in atom:
                continue
            chain = atom.split()[4].strip()
            if chain not in colourByChain:
                colourByChain[chain]='black'
    keysC = sorted(colourByChain.keys())
    print ('\n**************************************\n')
    print ("ARROW COLOUR KEY BY CHAIN")
    for k in keysC:
        print (k+': '+colourByChain[k])
    print ('\n**************************************\n')
    # Get the vectors
    # CHANGE HERE # may not be correct change, double check

    #Write VISUALISE
    # Control writing animations
    nsteps = 1
    if animate:
        nsteps = 50
    try:
        outpdbfilename = outdir + "/VISUALISE/" + structure + '_' + str(mode+1) + ".pdb"
        w = open(outpdbfilename, 'w')
        for i in range(0, nsteps):
            v_index = -1
            for atom in c_beta_atoms:
                if "ATOM" in atom:
                    v_index += 1
                    v = vectors[v_index]
                    vx = v[0]
                    vy = v[1]
                    vz = v[2]
                    x = round(float(atom[30:38].strip()) + (vx * i / 5), 3)
                    xspace = ' ' * (len(atom[30:38]) - len(str(x)))
                    y = round(float(atom[38:46].strip()) + (vy * i / 5), 3)
                    yspace = ' ' * (len(atom[38:46]) - len(str(y)))
                    z = round(float(atom[46:54].strip()) + (vz * i / 5), 3)
                    zspace = ' ' * (len(atom[46:54]) - len(str(z)))
                    atom = atom[0:30] + xspace + str(x) + yspace + str(y) +\
                                        zspace + str(z) + atom[54:]
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
                                    w.write("CONECT" + " " * (5 - len(atom1)) +\
                                            atom1 + " " * (5 - len(atom2)) + atom2 + "\n")
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
              '$start]]]\n    graphics $mol cylinder $start $middle radius '+str(tail)+'\n    graphics $mol cone $middle $end '
              'radius '+str(head)+'\n}\ndraw color black\n')
        else:
            arrows.append('proc vmd_draw_arrow {mol start end} {\n    set middle [vecadd $start [vecscale 0.9 [vecsub $end '
              '$start]]]\n    graphics $mol cylinder $start $middle radius '+str(tail)+'\n    graphics $mol cone $middle $end '
              'radius '+str(head)+'\n}\ndraw color '+colours[0]+'\n')
        recognised = False
        nextAtom=[]
        for n,atom in enumerate(c_beta_atoms,1):
            if "ATOM" in atom:
                chain = atom.split()[4].strip()
                v_index += 1
                if chain == chainSpecified or chainSpecified =='all':
                    recognised = True
                    v = vectors[v_index]
                    vx = v[0]
                    vy = v[1]
                    vz = v[2]
                    x1 = round(float(atom[30:38].strip()) + (vx * steps[0]), 3)
                    y1 = round(float(atom[38:46].strip()) + (vy * steps[0]), 3)
                    z1 = round(float(atom[46:54].strip()) + (vz * steps[0]), 3)
                    x2 = round(float(atom[30:38].strip()) + (vx * steps[1])*arrowLength, 3)
                    y2 = round(float(atom[38:46].strip()) + (vy * steps[1])*arrowLength, 3)
                    z2 = round(float(atom[46:54].strip()) + (vz * steps[1])*arrowLength, 3)
                    arrows.append('draw arrow {' + str(x1) + ' ' + str(y1) + ' ' + str(z1) + '} {' + str(x2) + ' ' + str(y2) + ' ' + str(z2) + '}\n')
                    if n in chainbreaks and chainSpecified == 'all':
                        nextAtom = c_beta_atoms[n:]
                        for a in nextAtom:
                            if "ATOM" in a:
                                chain = a.split()[4].strip()
                                break
                        arrows.append('draw color '+colourByChain[chain]+'\n')
        arrow_filename = outdir + "/VISUALISE/" + structure + '_ARROWS_' + str(mode+1) + ".txt"
        ngl_filename = "{0}/VISUALISE/{1}_NGL_ARROWS.js".format(outdir, mode+1)
        w = open(arrow_filename, 'w')
        w.writelines(arrows)
        w.close()
        args.filename = arrow_filename
        args.outfilename = ngl_filename
        visualisevector2ngl(args)
        if not recognised:
            print ('ERROR: Specified chain: "{}" not recognised'.format(chainSpecified))
    except IndexError:
        print ('ERROR: Vector file and pdb file are not compatible')
        sys.exit(1)
    except KeyError:
        print ('ERROR: MODE-TASK cannot colour all chains. Colours should match the number of chains')
        sys.exit(1)

def main(args):

    head = args.head
    tail = args.tail
    arrowLength = args.arrowLength
    coloursSpecified = args.colourByChain
    chainSpecified = args.chain
    aUnits = args.aUnits.split(',')
    atomT = args.atomType.upper()
    if atomT!='CA' and atomT!='CB':
        print('ERROR: Unrecognised atom type (Choose "CA" to select alpha carbons, or CB for beta carbons')
        sys.exit(1)
    try:
        pdb_file = open(args.pdb, 'r')
        pdb_lines = pdb_file.readlines()
        pdb_file.close()
    except IOError:
        print('ERROR: File {} not found'.format(args.pdb))
        sys.exit(1)

    # CHANGE THE FOLLOWING
#    mode = args.mode
    outdir = args.outdir
    vtMatrix = args.vtMatrix
    direction = args.direction
    if direction !=1 and direction !=-1:
        direction = 1
        print ('WARNING: Direction can only be 1 or -1: Default direction (+1) has been used')
    total_nmodes = get_nmodes(vtMatrix)
    if args.nmodes > total_nmodes:
        print("ERROR: Number of modes exceeds the actual number [{}]".format(total_nmodes))
        sys.exit(1)
    if args.nmodes < 1:
        print("ERROR: Choose a number grater than 0")
        sys.exit(1)
    for modenum in range(args.nmodes):
        vectorFile = getVector(vtMatrix, modenum, direction)
        visualise(pdb_lines, atomT, modenum, vectorFile, outdir, head, tail,
                  arrowLength, coloursSpecified, chainSpecified, aUnits,
                  animate=args.animate)

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
    parser.add_argument("--silent", action='store_true', default=False,
                        help="Turn off logging")
    parser.add_argument("--welcome", default="true",
                        help="Display welcome message (true/false)")
    parser.add_argument("--log-file", default=None,
                        help="Output log file (default: standard output)")
    parser.add_argument("--outdir", default="output",
                        help="Output directory")
    # custom arguments
    parser.add_argument("--pdb", help="Coarse grained PDB file")
    parser.add_argument("--nmodes", type=int, default=50,
                        help="Number of modes to select (default=50)")
    parser.add_argument("--direction", type=int, default = 1,
                        help="[int 1 or -1] Direction of overlap correction (Default = 1)")
    parser.add_argument("--vtMatrix",
                        help="File containing eigen vectors")
    parser.add_argument("--atomType", default='X',
                        help="Enter CA to select alpha carbons or CB to select beta carbons")
    # additional opitions
    parser.add_argument("--head", type = float, default=2.20,
                        help="Specify radius radius of arrow head. Default = 2.20")
    parser.add_argument("--tail", type = float, default=0.80,
                        help="Specify radius radius of arrow tail cylinder. Default = 0.80")
    parser.add_argument("--arrowLength",type=float, default=1.0,
                        help="Increase of decrease arrow length by a specified factor")
    parser.add_argument("--colourByChain", default='none', help="Comma separated list of colours for each chain.\nEg: If you have four chains in your complex then enter something like: red,blue,green,pink\nThis will colour the following:\nChain 1 as red\nChain 2 as blue\nChain 3 as green\nChain 4 as pink\nNote that these colours must match available colours in VMD\n")
    parser.add_argument("--chain", default='all',
                        help="Draw arrows for specified chain only")
    parser.add_argument("--aUnits", default='all',
                        help="Draw arrows for specified asymmetric units only. 1) Specify a single unit.\n  E.g --aUnit 1 \nOR\n 2) Provide a list of asymmetric units.\n  E.g --aUnits 1,4,5")
    parser.add_argument("--animate", action="store_true",
                        help="Add eigenvectors to structure multiple times to simulate motion")

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
        print ('No arguments provided. Use -h to view help')
