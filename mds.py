#!/usr/bin/env python
#filename: mds.py
import os
import sys
import argparse
import mdtraj as md
import numpy as np
from sklearn.metrics import euclidean_distances
from sklearn.manifold import MDS
from lib.write_plot import write_plots, write_pcs, write_fig
from lib.traj_info import trajectory_info, get_internal_coordinates, get_trajectory, get_cosine, print_kmo
from lib.utils import welcome_msg

def main():
	
	return;

#==============================================================================#
#											MDS MD
#
#			This programe performs the MDS on a MD trajectory
#
# 								Author : Bilal Nizami
# 						  	 Rhodes University, 2017
#==============================================================================#

##===============================================================================
##								 Welcome message
##===============================================================================
title='MDS MD'
welcome_msg(title, "Bilal Nizami")

#==============================================================================
#                            Setting the options
#==============================================================================

def get_options():
	parser = argparse.ArgumentParser(usage='%(prog)s -t <MD trajectory> -p <topology file>')
	
	parser.add_argument("-t", "--trj", dest="trj",				help="file name of the MD trajectory")
	parser.add_argument("-p", "--top", dest="topology",				help="topology file")
	parser.add_argument("-out", "--out", dest="out_dir", help="Name of the output directory. Default is out")
	parser.add_argument("-mt", "--mds_type", dest="mds_type",				help="Type of MDS. Options are nm=non-metric, metric=metric")
	parser.add_argument("-dt", "--dissimilarity_type",  dest="dissimilarity_type",				help="Type of dissimilarity matrix to use. euc = Euclidean distance between internal coordinates, rmsd= pairwise RMSD. Default is rmsd")
	parser.add_argument("-ag", "--ag", dest="atm_grp", help="group of atom for MDS. Default is C alpha atoms. Other options are :"				  "all= all atoms, backbone = backbone atoms, CA= C alpha atoms, protein= protein's atoms")	
	parser.add_argument("-ct", "--coordinate_type",  dest="coordinate_type",				help="Internal coordinates type. Default is pairwise distance")
	parser.add_argument("-ai", "--atom_indices",  dest="atom_indices",				help="group of atom for pairwise distance. Default is C alpha atoms. Other options are :"				  "all= all atoms, backbone = backbone atoms, alpha= C alpha atoms, heavy= all non hydrogen atoms, minimal=CA,CB,C,N,O atoms")

	args = parser.parse_args()
	
	if args.out_dir == None:
		out=args.trj
		args.out_dir=out
		
	if args.trj is None:
		print 'ERROR: Missing trajectory argument.... :(  \nPlease see the help by running \n\nsystem_setup.py -h\n\n '
		parser.print_help()
		sys.exit(1)

	if not os.path.exists(args.trj ):
		print('\nERROR: {0} not found....:(  Please check the path or filename\n' .format(args.trj ))
		#parser.print_help()
		sys.exit(1)
		
	if not os.path.exists(args.topology):
		print('\nERROR: {0} not found....:(  Please check the path or filename\n' .format(args.topology ))
		#parser.print_help()
		sys.exit(1)
		
	if args.topology is None:
		print 'ERROR: Missing toplogy argument.... :(  \nPlease see the help by running \n mds.py -h\n\n '
		parser.print_help()
		sys.exit(1)
		
	if args.dissimilarity_type not in  ('rmsd', 'euc', None):
		print 'ERROR: no such option as', args.dissimilarity_type, 'for flag -dt \nPlease see the help by running \n mds.py -h\n\n '
		sys.exit(1)
		
	if args.mds_type == None:
		print 'No MDS type given... performing metric MDS...'
		args.mds_type='metric'
	if args.mds_type not in  ('nm', 'metric', None):
		print 'ERROR: no such option as', args.mds_type, 'for flag -mt \nPlease see the usage\n\n '
		sys.exit(1)
		
	if args.coordinate_type not in  ('distance', 'phi', 'psi', 'angle', None):
		print 'ERROR: no such option as', args.coordinate_type, 'for flag -ct \nPlease see the usage\n\n '
		parser.print_help()
		sys.exit(1)
	
	if args.atm_grp == None:
		print 'No atom selected. MDS will be performed on C alpha atoms '
		args.atm_grp = 'CA'  # set to default C-alpha atoms
		
	if args.atm_grp not in  ('all', 'CA', 'backbone', 'protein'):
		print 'ERROR: no such option as', args.atm_grp, 'for flag -at \nPlease see the usage\n\n '
		sys.exit(1)
	if 	args.dissimilarity_type == 'euc':
		if args.atom_indices == None:
			print 'No atom selected for pairwise distance. pairwise distance of C alpha atoms will be used'
			args.atom_indices='alpha'
	
	if args.atom_indices not in  ('all', 'alpha', 'backbone', 'minimal', 'heavy', None):
		print 'ERROR: no such option as', args.atom_indices, 'for flag -ai \nPlease see the usage\n\n '
		sys.exit(1)
	if args.dissimilarity_type == None or args.dissimilarity_type == 'rmsd':
		if args.atom_indices != None:
			print '\nWARNING: -ai', args.atom_indices, '  ,is meaningless with -dt set to rmsd \n'
	
	return args;

args=get_options()

#=======================================
# assign the passed arguments and read the trajectory
#=======================================

traj = args.trj
topology = args.topology
#pca_traj = md.load(traj, top=topology)

print 'Reading trajectory ', args.trj, '...' 

try:
	pca_traj = md.load(traj, top=topology)
except:
	raise IOError('Could not open trajectory {0} for reading. \n' .format(trj))
top = pca_traj.topology
atm_name=args.atm_grp
sele_grp = get_trajectory(atm_name, top)
atom_indices=args.atom_indices

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
## =============================
# print trajectory info
#===================================
trajectory_info(pca_traj, traj, atm_name, sele_grp)

# print KMO 
print_kmo(pca_traj, traj, atm_name, sele_grp)

def get_pair_rmsd(pca_traj, sele_grp):
	'pair wise RMSD over all the frames, return a square matrix of pairwise rmsd'
	pair_rmsd=np.empty((pca_traj.n_frames, pca_traj.n_frames))
	for i in range(pca_traj.n_frames):
		pair_rmsd[i]=md.rmsd(pca_traj, pca_traj, i, atom_indices=sele_grp)
	pair_rmsd=(pair_rmsd+pair_rmsd.transpose())/2  ## due to precision level matrix might not evaluate as symmetric, hence to make it symmetric
	return pair_rmsd;


#============================================
#
#  Multidimensional scaling
#
#=============================================

def mds(input, type):
	'metric and nonmetric Multidimensional scaling'
	seed = np.random.RandomState(seed=1)
	#np.savetxt('mds_input.txt', input) ## testing value error
	if type == 'nm':
		nmds=MDS(n_components=100, max_iter=3000, metric=False, random_state=seed, dissimilarity="precomputed")
		print "Performing non-metric MDS.."
		npos=nmds.fit_transform(input)
		
		# write PC plots
		write_plots('nmds_projection', npos, out_dir)
		title='nMDS Projection'
		write_fig('nmds_projection', npos, out_dir, title)
		# cosine content 
		
		pc1_cos=get_cosine(npos, 0)
		print 'cosine content of first PC=',pc1_cos
		pc2_cos=get_cosine(npos, 1)
		print 'cosine content of second PC=', pc2_cos
		pc3_cos=get_cosine(npos, 2)
		print 'cosine content of 3rd PC=',pc3_cos
		pc4_cos=get_cosine(npos, 3)
		print 'cosine content of 4th PC=', pc4_cos
			
	elif type == 'metric':
		mmds = MDS(n_components=100, max_iter=3000, random_state=seed, dissimilarity="precomputed")
		print "Performing metric MDS.."
		mpos = mmds.fit_transform(input)
		
		# write PC plots
		write_plots('mmds_projection', mpos, out_dir)
		title='mMDS Projection'
		write_fig('mmds_projection', mpos, out_dir, title)
		# cosine content 
		
		pc1_cos=get_cosine(mpos, 0)
		print 'cosine content of first PC=',pc1_cos
		pc2_cos=get_cosine(mpos, 1)
		print 'cosine content of second PC=', pc2_cos
		pc3_cos=get_cosine(mpos, 2)
		print 'cosine content of 3rd PC=',pc3_cos
		pc4_cos=get_cosine(mpos, 3)
		print 'cosine content of 4th PC=', pc4_cos
	
	else:
		print 'ERROR: Please check -mt flag options by running mds.py -h'
	
	return;

type=args.mds_type

if args.dissimilarity_type == 'rmsd' or args.dissimilarity_type == None:
	print 'using pairwise RMSD...'
	pair_rmsd=get_pair_rmsd(pca_traj, sele_grp)
	mds(pair_rmsd, type)
	print 'FINISHED!'
if args.dissimilarity_type == 'euc':
	if args.coordinate_type == None:
		args.coordinate_type = "distance"
		print "Using pairwise distance by default"
	print 'using Euclidean space of', args.coordinate_type
	int_cord=get_internal_coordinates(top, args.coordinate_type, pca_traj, atom_indices)
	similarities = euclidean_distances(int_cord)
	mds(similarities, type)
	print 'FINISHED!'



if __name__=="__main__":
	main()