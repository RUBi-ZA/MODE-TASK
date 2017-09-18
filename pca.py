#!/usr/bin/env python
#filename: pca.py
import os
import sys
import time
import re
import math
import argparse
from time import sleep, gmtime, strftime
from datetime import datetime
import mdtraj as md
import numpy as np
from sklearn.decomposition import PCA, KernelPCA, IncrementalPCA
from sklearn.metrics import euclidean_distances
from sklearn import preprocessing
from lib.write_plot import write_plots, write_pcs, write_fig
from lib.traj_info import trajectory_info, get_cosine, get_kmo, print_kmo, get_rmsf
from lib.utils import welcome_msg
import scipy.integrate

def main():
	
	return;


#==============================================================================#
#											PCA MD 
#
#			This programe performs the PCA on a MD trajectory
#
# 								Author : Bilal Nizami
# 						  	 Rhodes University, 2017
#==============================================================================#

##===============================================================================
##								 Welcome message
##===============================================================================
title='PCA MD'
welcome_msg(title, "Bilal Nizami")


#==============================================================================
#                            Setting the options
#==============================================================================

def set_option():
	parser = argparse.ArgumentParser( usage='%(prog)s -t <MD trajectory> -p <topology file>')
	#"Usage: pca.py -t <MD trajectory> -p <topology file>  -a <atom group >"
	
	parser.add_argument("-t", "--trj", dest="trj", help="file name of the MD trajectory", action="store")
	parser.add_argument("-p", "--top", dest="topology", help="topology file")     
	parser.add_argument("-out", "--out", dest="out_dir", help="Name of the output directory. Default is out")
	parser.add_argument("-ag", "--ag", dest="atm_grp", help="group of atom for PCA. Default is C alpha atoms. Other options are :"				  "all= all atoms, backbone = backbone atoms, CA= C alpha atoms, protein= protein's atoms")	
	parser.add_argument("-r", "--ref", dest="reference", help="reference structure for RMSD") 
	parser.add_argument("-pt", "--pca_type", dest="pca_type", help="PCA method. Default is svd (Single Value Decomposition) PCA. Options are:\
					evd, kpca, svd, ipca. If svd is selected, additional arguments can be passed by flag -st. If KernelPCA is selected kernel type can also be defined by flag -kt") 	
	parser.add_argument("-nc", "--comp", type=int, dest="comp", help="Number of components to keep in a PCA object. If not set, by default all the components will be kept.")	
	parser.add_argument("-kt", "--kernel_type", dest="kernel_type", help="Type of kernel for KernelPCA. default is linear. Options are :"
					"linear, poly, rbf, sigmoid, cosine, precomputed") 
	parser.add_argument("-st", "--svd_solver", dest="svd_solver", help="Type of svd_solver for SVD (Single Value Decomposition) PCA. Default is auto. Options are :"				  "auto, full, arpack, randomized") 
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
	if args.pca_type not in  ('svd', 'evd', 'kpca', 'ipca', None):
		print 'ERROR: no such option as', args.pca_type, 'for flag -pt \nPlease see the help by running \n pca.py -h..\n\n '
		sys.exit(1)
	if args.kernel_type not in  ('linear', 'poly', 'sigmoid', 'cosine', 'precomputed', 'rbf',None):
		print 'ERROR: no such option as', args.kernel_type, 'for flag -kt \nPlease see the help by running \n pca.py -h..\n\n '
		sys.exit(1)
		
	if args.kernel_type != None and args.pca_type != 'kpca':
		print 'WARNING: -kt', args.kernel_type, 'is meaningless with -pt', args.pca_type, '. Flag -kt is being ignored!'
		
	if args.svd_solver not in  ('auto', 'full', 'arpack', 'randomized', None):
		print 'ERROR: no such option as', args.svd_solver, 'for flag -st \nPlease see the help by running \n pca.py -h.\n\n'
		sys.exit(1)
	
	if args.svd_solver != None and args.pca_type != 'svd':
		print 'WARNING: -st', args.svd_solver, 'is meaningless with -pt', args.pca_type, '. Flag -st is being ignored!'
	return args

args = set_option()
atm_name = args.atm_grp
#=======================================
# assign the passed arguments and read the trajectory 
#=======================================
traj = args.trj
topology = args.topology
ref = args.reference
ptype=args.pca_type

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

comp = args.comp
n_eivec=5
## read the reference structure
if ref:
	try:
		ref = md.load(args.reference)
	except:
			raise IOError('Could not open reference structure {0} for reading. \n' .format(args.reference))

print 'Reading trajectory ', args.trj, '...' 
try:
	pca_traj = md.load(traj, top=topology)
except:
	raise IOError('Could not open trajectory {0} for reading. \n' .format(trj))
	
top = pca_traj.topology


#==============================================
#
# Setting the default options
#
#===============================================

if args.atm_grp == None:
	print 'No atom has been selected. PCA will be performed on C alpha atoms '
	atm_name = 'CA'  # set to default C-alpha atoms
if args.reference == None:
	print "No reference structure given, RMSD will be computed to the first frame in the trajectory"
	ref = pca_traj # set reference to current trajectory
if args.pca_type == None:
	ptype = 'svd'

if args.svd_solver == None:
	svd='auto'
	
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

# Print trajectory information
trajectory_info(pca_traj, traj, atm_name, sele_grp)

# print KMO 
#sele_grp1=top.select("backbone") ## select protein for KMO
print_kmo(pca_traj, traj, atm_name, sele_grp)

#===============================================================
#
#  RMSD in reference with first frame
#
#===============================================================

def get_rmsd():
	rmsd = md.rmsd(pca_traj, ref, 0, atom_indices=sele_grp)
	print "RMSD written to rmsd.agr \n "
	## write the RMSD file
	rmsd_fname=out_dir+'/rmsd.agr'
	np.savetxt(rmsd_fname, rmsd)
	rf = open(rmsd_fname, 'r')
	rf_cont = rf.read()
	rf.close()
	
	my_time = strftime("%Y-%m-%d  %a  %H:%M:%S", gmtime())
	title = '\tcreated by pca.py\t'
	legends = '@    title "RMSD"\n\
	@    xaxis  label "Time (ps)"\n\
	@    yaxis  label "RMSD (\cE\C) "\n\
	@	TYPE xy\n'
	
	pf = open(rmsd_fname, 'w')
	pf.write('#'+title+'\ton\t'+my_time+'\n'+legends+'\n'+rf_cont)
	pf.close()
	return;
	
get_rmsd()

#===============================================================
#
# PCA using sci-kit learn library
#===============================================================

def svd_pca(svd):
	'single value decomposition based PCA'
	pca_traj.superpose(pca_traj, 0, atom_indices=sele_grp) 			# Superpose each conformation in the trajectory upon first frame
	sele_trj = pca_traj.xyz[:,sele_grp,:]												# select coordinates of selected atom groups
	sele_traj_reshaped = sele_trj.reshape(pca_traj.n_frames, len(sele_grp) * 3)
	sele_traj_reshaped = sele_traj_reshaped.astype(float) ## to avoid numpy Conversion Error during scaling
	sele_traj_reshaped_scaled = preprocessing.scale(sele_traj_reshaped, axis=0, with_std=False) # center to the mean
	
	pca_sele_traj = PCA(n_components=comp)
	pca_sele_traj.fit(sele_traj_reshaped_scaled)
	pca_sele_traj_reduced = pca_sele_traj.transform(sele_traj_reshaped_scaled)
	
	print "Trace of the covariance matrix is: ", np.trace(pca_sele_traj.get_covariance())
	
	# write the plots 
	
	write_plots('pca_projection', pca_sele_traj_reduced, out_dir)
	title='PCA Projection'
	write_fig('pca_projection', pca_sele_traj_reduced, out_dir, title)

	#write the pcs variance
	write_pcs('pca_variance', pca_sele_traj, out_dir)
	
	pc1_cos=get_cosine(pca_sele_traj_reduced, 0)
	print 'cosine content of first PC=',pc1_cos
	pc2_cos=get_cosine(pca_sele_traj_reduced, 1)
	print 'cosine content of second PC=', pc2_cos
	pc3_cos=get_cosine(pca_sele_traj_reduced, 2)
	print 'cosine content of 3rd PC=',pc3_cos
	pc4_cos=get_cosine(pca_sele_traj_reduced, 3)
	print 'cosine content of 4th PC=', pc4_cos
	
	return;


#==============================================================
#
# Kernel PCA
#
# ==============================================================
def my_kernelPCA(kernel):
	pca_traj.superpose(pca_traj, 0, atom_indices=sele_grp) 			# Superpose each conformation in the trajectory upon first frame
	sele_trj = pca_traj.xyz[:,sele_grp,:]												# select coordinates of selected atom groups
	sele_traj_reshaped = sele_trj.reshape(pca_traj.n_frames, len(sele_grp) * 3)
	sele_traj_reshaped = sele_traj_reshaped.astype(float) ## to avoid numpy Conversion Error during scaling
	sele_traj_reshaped_scaled = preprocessing.scale(sele_traj_reshaped, axis=0, with_std=False) # center to the mean

	kpca = KernelPCA(kernel = kernel, fit_inverse_transform=True, gamma=10)
	kpca.fit(sele_traj_reshaped_scaled)
	#print "Trace of the covariance matrix is: ", np.trace(kpca.get_covariance())
	kpca_reduced = kpca.transform(sele_traj_reshaped_scaled)
	
	#write plots
	write_plots('kpca_projection', kpca_reduced, out_dir)
	title='kPCA Projection'
	write_fig('kpca_projection', kpca_reduced, out_dir, title)

	#write variance
	kpca_variance_fname=out_dir+'/kpca_variance'
	np.savetxt(kpca_variance_fname, kpca.lambdas_)
	
	pc1_cos=get_cosine(kpca_reduced, 0)
	print 'cosine content of first PC=',pc1_cos
	pc2_cos=get_cosine(kpca_reduced, 1)
	print 'cosine content of second PC=', pc2_cos
	pc3_cos=get_cosine(kpca_reduced, 2)
	print 'cosine content of 3rd PC=',pc3_cos
	pc4_cos=get_cosine(kpca_reduced, 3)
	print 'cosine content of 4th PC=', pc4_cos
	return;


#=============================================================
#
# Incremental PCA
#
#=============================================================

def incremental_pca():
	' normal PCA is very memory intesive. It can be problemetic for large dataset, \
	since dataset is stored in memory. Incremental principal component analysis (IPCA) is \
	typically used for such cases. '
	
	pca_traj.superpose(pca_traj, 0, atom_indices=sele_grp) 			# Superpose each conformation in the trajectory upon first frame
	sele_trj = pca_traj.xyz[:,sele_grp,:]												# select coordinates of selected atom groups
	sele_traj_reshaped = sele_trj.reshape(pca_traj.n_frames, len(sele_grp) * 3)
	sele_traj_reshaped = sele_traj_reshaped.astype(float) ## to avoid numpy Conversion Error during scaling
	sele_traj_reshaped_scaled = preprocessing.scale(sele_traj_reshaped, axis=0, with_std=False) # center to the mean

	ipca = IncrementalPCA()
	ipca = ipca.fit(sele_traj_reshaped_scaled)
	ipca_reduced=ipca.transform(sele_traj_reshaped_scaled)
	
	#write plots
	write_plots('ipca_projection', ipca_reduced, out_dir)
	title='iPCA Projection'
	write_fig('ipca_projection', ipca_reduced, out_dir, title)
	#write variance
	#np.savetxt('ipca_variance', kpca.lambdas_)
	pc1_cos=get_cosine(ipca_reduced, 0)
	print 'cosine content of first PC=',pc1_cos
	pc2_cos=get_cosine(ipca_reduced, 1)
	print 'cosine content of second PC=', pc2_cos
	pc3_cos=get_cosine(ipca_reduced, 2)
	print 'cosine content of 3rd PC=',pc3_cos
	pc4_cos=get_cosine(ipca_reduced, 3)
	print 'cosine content of 4th PC=', pc4_cos

	return;

#===============================================================
#
#  Eigenvalue decomposition based PCA
#
#=================================================================
def my_pca():
	'eigenvales decomposition PCA'
	pca_traj.superpose(pca_traj, 0, atom_indices=sele_grp) 			# Superpose each conformation in the trajectory upon first frame
	sele_trj = pca_traj.xyz[:,sele_grp,:]												# select coordinates of selected atom groups
	sele_traj_reshaped = sele_trj.reshape(pca_traj.n_frames, len(sele_grp)* 3)
	#arr1=sele_traj_reshaped
	sele_traj_reshaped = sele_traj_reshaped.astype(float) ## to avoid numpy Conversion Error during scaling
	sele_traj_reshaped_scaled = preprocessing.scale(sele_traj_reshaped, axis=0, with_std=False) # center to the mean
	arr = sele_traj_reshaped_scaled
	
	#===============================================
	# covariance matrix 
	cov_mat = np.corrcoef(arr, rowvar=False)
	trj_eval, trj_evec=np.linalg.eig(cov_mat)
	
	print "Trace of cov matrix is ",  np.trace(cov_mat)
	
	#=============================
	# sanity check of calculated eigenvector and eigen values 
	# it must be cov matrix * eigen vector = eigen vector * eigen value
	
	for i in range(len(trj_eval)):
		eigv = trj_evec.real[:,i].reshape(1,len(trj_evec[:,0]),).T
		np.testing.assert_array_almost_equal(cov_mat.dot(eigv), trj_eval[i]*eigv, decimal=3, err_msg='', verbose=True)

#=============================================
	# sort the eigenvalues and eigenvector
	sort_idx = trj_eval.argsort()[::-1]
	trj_eval = trj_eval[sort_idx]
	trj_evec = trj_evec[sort_idx]

	tot_var = np.sum(trj_eval.real)
	variation = []
	cum = []
	j = 0
	eigv = []
	n_comp=100
	pca = trj_evec.real[:,0:n_comp]    ## keep first 100 eigenvectors
	for i in trj_eval.real[0:n_comp]:
		eigv.append(i)
		variation.append((i/tot_var)*100)
		j +=1
	
	# write PC plot 
	pca_variance_fname=out_dir+'/pca_variance.agr'
	np.savetxt(pca_variance_fname, variation)
	ef = open(pca_variance_fname, 'r')
	ef_cont = ef.read()
	ef.close()
	title = '\tcreated by pca.py\t'
	my_time = strftime("%Y-%m-%d  %a  %H:%M:%S", gmtime())
	legends = '@    title "explained_variance of PCs"\n\
	@    xaxis  label "PCs"\n\
	@    yaxis  label "% Variance"\n\
	@	TYPE xy\n\
	@    s0 symbol 1\n\
	@    s0 symbol size 0.250000\n\
	@    s0 symbol color 1\n\
	@    s0 symbol pattern 1\n\
	@    s0 symbol fill color 1\n\
	@    s0 symbol fill pattern 1\n\
	@    s0 symbol linewidth 1.0\n\
	@    s0 symbol linestyle 1\n\
	@    s0 symbol char 25\n\
	@	s0 symbol fill color 2\n\
	@	s0 symbol color 2\n\
	@    s0 symbol char font 0\n\
	@    s0 symbol skip 0\n'
	
	ef = open(pca_variance_fname, 'w')
	ef.write('#'+title+'\ton\t'+my_time+'\n'+legends+'\n'+ef_cont+'\n')
	ef.close()
	#========================================================
	# transform the input data into choosen pc
	arr_transformed = arr.dot(pca)
	#arr_transformed = pca.T.dot(arr.T)
	print arr_transformed.shape
	write_plots('pca_projection', arr_transformed, out_dir)
	title='PCA Projection'
	write_fig('pca_projection', arr_transformed, out_dir, title)
	
	## RMSF 
	get_rmsf(pca_traj, sele_grp, trj_eval, out_dir)
	
	pc1_cos=get_cosine(arr_transformed, 0)
	print 'cosine content of first PC=',pc1_cos
	pc2_cos=get_cosine(arr_transformed, 1)
	print 'cosine content of second PC=', pc2_cos
	pc3_cos=get_cosine(arr_transformed, 2)
	print 'cosine content of 3rd PC=',pc3_cos
	pc4_cos=get_cosine(arr_transformed, 3)
	print 'cosine content of 4th PC=', pc4_cos
	
	return;



if ptype == 'kpca':
	kernel = ''
	kernel = args.kernel_type
	if args.kernel_type:
		print "Performing Kernel PCA with", kernel, 'kernel'
		my_kernelPCA(kernel)
		print "\nFINISHED. !"
	else:
		print "Performing Kernel PCA with default linear kernel"
		my_kernelPCA('linear')
		print "\nFINISHED. !"

if ptype == 'svd':
	svd=''
	svd = args.svd_solver
	if svd:
		print "Performing SVD (Single Value Decomposition) PCA with ",svd,"svd_solver"
		svd_pca(svd)
		print "\nFINISHED. !"
	else:
		print "Performing SVD (Single Value Decomposition) PCA with 'auto' svd_solver"
		svd_pca(svd)
		print "\nFINISHED. !"
if ptype == 'ipca':
	print "Performing Incremental_pca (IPCA)"
	incremental_pca()
	print "\nFINISHED. !"

if ptype == 'evd':
	my_pca()
	print "\nFINISHED. !"

if __name__=="__main__":
	main()