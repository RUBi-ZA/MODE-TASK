#!/usr/bin/env python
#filename: internal_pca.py

import os
import sys
import argparse
import mdtraj as md
import numpy as np
from sklearn.decomposition import PCA, KernelPCA, IncrementalPCA
from sklearn.metrics import euclidean_distances
from itertools import combinations
from lib.write_plot import write_plots, write_pcs, write_fig
from lib.traj_info import trajectory_info, get_cosine, print_kmo
from lib.utils import welcome_msg

def main():
	
	return;
#==============================================================================#
#											Internal PCA MD 
#
#			This programe performs the PCA on internal coordinates of a MD trajectory
#
# 								Author : Bilal Nizami
# 						  	 Rhodes University, 2017
#==============================================================================#

##===============================================================================
##								 Welcome message
##===============================================================================

title='Internal PCA MD'
welcome_msg(title, "Bilal Nizami")

def set_option():
	parser = argparse.ArgumentParser( usage='%(prog)s -t <MD trajectory> -p <topology file>')
	#"Usage: pca.py -t <MD trajectory> -p <topology file>  -a <atom group >"
	
	parser.add_argument("-t", "--trj", dest="trj", help="file name of the MD trajectory", action="store")
	parser.add_argument("-p", "--top", dest="topology", help="topology file")   
	parser.add_argument("-out", "--out", dest="out_dir", help="Name of the output directory. Default is out")	
	parser.add_argument("-nc", "--comp", type=int, dest="comp", help="Number of components to keep in a PCA object. If not set, by default all the components will be kept.")	
	parser.add_argument("-ag", "--ag", dest="atm_grp", help="group of atom for PCA. Default is C alpha atoms. Other options are :"				  "all= all atoms, backbone = backbone atoms, CA= C alpha atoms, protein= protein's atoms")	
	parser.add_argument("-ct", "--ref", dest="coordinate_type", help="internal coordinate type. Options are: distance, angles, phi and, psi") 
	args = parser.parse_args()	
	
	
	#====================================================================
	# if no arguments are passed
	#====================================================================
	if args.out_dir == None:
		out=args.trj
		args.out_dir=out
		
	if args.trj is None: 
		print 'ERROR: Missing trajectory argument.... :(  \nPlease see the help by running \n\nsystem_setup.py -h\n\n '
		parser.print_help()
		sys.exit(1)
	
	if args.topology is None:
		print 'ERROR: Missing topology.... :( \nPlease see the help by running \n\nsystem_setup.py -h\n\n '
		parser.print_help()
		sys.exit(1)
	
	if not os.path.exists(args.trj ):
		print('\nERROR: {0} not found....:(  Please check the path\n' .format(args.trj ))
		parser.print_help()
		sys.exit(1)
	
	if not os.path.exists(args.topology):
		print('\nERROR: {0} not found....:(  Please check the path\n' .format(args.topology ))
		parser.print_help()
		sys.exit(1)
	if args.coordinate_type is None:
		print "No arguments given for -ct...using distance as internal coordinate\n"
		args.coordinate_type='distance'
	return args
	
args = set_option()
atm_name = args.atm_grp
#====================================================================
# if no arguments are passed
#====================================================================
if args.trj is None: 
	print 'Missing trajectory arguments :(\nPlease see the help by running \n\nsystem_setup.py -h\n\n '
	parser.print_help()
	sys.exit(1)

if args.topology is None:
	print 'Missing topology !!\nPlease see the help by running \n\nsystem_setup.py -h\n\n '
	parser.print_help()
	sys.exit(1)


#=======================================
# assign the passed arguments and read the trajectory 
#=======================================

traj = args.trj
topology = args.topology
#pca_traj = md.load(traj, top=topology)
comp = args.comp
print 'Reading trajectory ', args.trj, '...' 

try:
	pca_traj = md.load(traj, top=topology)
except:
	raise IOError('Could not open trajectory {0} for reading. \n' .format(trj))
top = pca_traj.topology

# take the input trj name for output directory
out_dir=args.out_dir
out_dir=out_dir.split('/')
out_dir=out_dir[-1]
out_dir='out_'+out_dir

if not os.path.exists(out_dir):
    os.makedirs(out_dir)
#else:
#	print out_dir, 'already exist. Can not overwrite the output directory!\n'
#	sys.exit(1)
print 'Results will be written in ', out_dir

#==============================================
#
# Setting the default options
#
#===============================================

if args.atm_grp == None:
	print 'No atom has been selected. PCA will be performed on C alpha atoms '
	atm_name = 'CA'  # set to default C-alpha atoms


#==========================================================================
#
#		selecting the atoms 
#
# User passes the arguements to select the subset of atoms for PCA
#===========================================================================


def get_trajectory():
	'get the part of system for PCA based on users input of atom group'
	if atm_name == 'CA':
		sele_grp=top.select("name CA")	

	if atm_name == 'backbone':
		sele_grp=top.select("backbone")
		
	if atm_name == 'all':
		sele_grp=top.select("all")
	return sele_grp;
	
	if atm_name == 'protein':
		sele_grp=top.select("protein")
	return sele_grp;

sele_grp = get_trajectory()

# print trajectory informations	
trajectory_info(pca_traj, traj, atm_name, sele_grp)

# print KMO 
print_kmo(pca_traj, traj, atm_name, sele_grp)


#===========================================================
#
# Internal coordinate type
#
#===========================================================
def get_internal_coordinates():
	'get the different types of internal coordinates as per user selections'
	calpha_idx=top.select_atom_indices('alpha')
	if args.coordinate_type == 'distance':
		print "Pair wise atomic distance selected\n "
		atom_pairs = list(combinations(calpha_idx, 2)) # all unique pairs of elements 
		pairwise_distances = md.geometry.compute_distances(pca_traj, atom_pairs)
		int_cord=pairwise_distances
		#print int_cord.shape
	if args.coordinate_type == 'phi':
		print  "phi torsions  selected\n"
		atom_pairs = list(combinations(calpha_idx, 3)) 
		angle=md.compute_phi(pca_traj)
		
		int_cord=angle[1] ## apparently compute_phi returns tupple of atoms indices and phi angles, index 1 has phi angles 
		#print np.array(angle[1]).shape
		#print int_cord[0]
	
	if args.coordinate_type == 'psi':
		print "psi torsions  selected\n"
		atom_pairs = list(combinations(calpha_idx, 3)) 
		angle=md.compute_psi(pca_traj)
		
		int_cord=angle[1] ## apparently compute_psi returns tupple of atoms indices and psi angles, index 1 has psi angles 
		#print np.array(angle[1]).shape
		
	if args.coordinate_type == 'angle':
		print "1-3 angle selected between N,CA and C"
		nrow=len(top.select("name CA")) # to get the number of amino acid ignoring ligand etc. 
		ncol=3
		# make a matrix of N,CA, C index, each row index making bond
		B = np.ones((nrow, ncol))
		B[:,0]= top.select('backbone and name N')
		B[:,1]= top.select('backbone and name CA')
		B[:,2]= top.select('backbone and name C')
		# compute angle between N,CA, C
		angle=md.compute_angles(pca_traj, B)
		int_cord=angle
	return int_cord;

#===========================================================
#
#  Internal Distance Coordinate Based PCA
#
#===========================================================
def distance_pca(int_cord1):
	'Internal Coordinate Based PCA'
	
	pca = PCA(n_components=comp)
	dpca = pca.fit(int_cord1)
	dpca_reduced=dpca.transform(int_cord1)
	
	write_plots('dpca_projection', dpca_reduced, out_dir)
	write_pcs('dpca_pcs', dpca, out_dir)
	title='internal coordinate PCA Projection'
	write_fig('dpca_projection', dpca_reduced, out_dir, title)
	
	pc1_cos=get_cosine(dpca_reduced, 0)
	print 'cosine content of first PC=',pc1_cos
	pc2_cos=get_cosine(dpca_reduced, 1)
	print 'cosine content of second PC=', pc2_cos
	pc3_cos=get_cosine(dpca_reduced, 2)
	print 'cosine content of 3rd PC=',pc3_cos
	pc4_cos=get_cosine(dpca_reduced, 3)
	print 'cosine content of 4th PC=', pc4_cos
	return;

int_cord=get_internal_coordinates()
distance_pca(int_cord)

if __name__=="__main__":
	main()