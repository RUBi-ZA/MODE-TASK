#!/usr/bin/python
#filename: pca.py
import numpy as np
from time import sleep, gmtime, strftime

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt


#==============================================================================#
#											
#			This programe is the part of PCA MD. It writes the PCA plots in xmgrace formatted .agr file 
#
# 								Author : Bilal Nizami
# 						  	 Rhodes University, 2017
#==============================================================================#

## write plots
def write_plots(file_name, pca, out_dir):
	'function to write pca plots. takes name of the file to write and pca object name'
	my_time = strftime("%Y-%m-%d  %a  %H:%M:%S", gmtime())
	title = '\tcreated by pca.py\t'
	legends12 = '@    title "Projection of PC"\n\
	@    xaxis  label "PC1"\n\
	@    yaxis  label "PC2"\n'
	legends13 = '@    title "Projection of PC"\n\
	@    xaxis  label "PC1"\n\
	@    yaxis  label "PC3"\n'
	legends23 = '@    title "Projection of PC"\n\
	@    xaxis  label "PC2"\n\
	@    yaxis  label "PC3"\n'
	other='@	TYPE xy\n\
	@    s0 line type 0\n\
	@    s0 line linestyle 1\n\
	@    s0 line linewidth 1.0\n\
	@    s0 line color 1\n\
	@    s0 line pattern 1\n\
	@    s0 baseline type 0\n\
	@    s0 baseline off\n\
	@    s0 dropline off\n\
	@    s0 symbol 1\n\
	@    s0 symbol size 0.250000\n\
	@    s0 symbol color 1\n\
	@    s0 symbol pattern 1\n\
	@    s0 symbol fill color 1\n\
	@    s0 symbol fill pattern 1\n\
	@    s0 symbol linewidth 1.0\n\
	@    s0 symbol linestyle 1\n\
	@    s0 symbol char 25\n\
	@    s0 symbol char font 0\n\
	@    s0 symbol skip 0\n'
	
	fname = ''
	
	#=======================================
	# xmgrace formated plot
	#================================
	# plot 1 and 2 PC
	
	pca1=pca[:,[0,1]]
	fname = out_dir+'/'+file_name+'1_2'+'.agr'
	np.savetxt(fname, pca1)
	pf = open(fname, 'r')
	pf_cont = pf.read()
	pf.close()
	pf = open(fname, 'w')
	pf.write('#'+title+'\ton\t'+my_time+'\n'+legends12+other+'\n'+pf_cont)
	pf.close()
	
	# plot 1 and 3 PC
	pca1=pca[:,[0,2]]
	fname = out_dir+'/'+file_name+'1_3'+'.agr'
	np.savetxt(fname, pca1)
	pf = open(fname, 'r')
	pf_cont = pf.read()
	pf.close()
	pf = open(fname, 'w')
	pf.write('#'+title+'\ton\t'+my_time+'\n'+legends13+other+'\n'+pf_cont)
	pf.close()
	
	# plot 2 and 3 PC
	pca1=pca[:,[1,2]]
	fname = out_dir+'/'+file_name+'2_3'+'.agr'
	np.savetxt(fname, pca1)
	pf = open(fname, 'r')
	pf_cont = pf.read()
	pf.close()
	pf = open(fname, 'w')
	pf.write('#'+title+'\ton\t'+my_time+'\n'+legends23+other+'\n'+pf_cont)
	pf.close()
	
	return;

def write_fig(file_name, pca, out_dir, title):
	#========================================
	# coloured png using matplotlib
	#=================================================
	## # plot 1 and 2 PC
	fname = out_dir+'/'+file_name+'1_2'+'.png'
	fig=plt.figure()
	col=range(1,len(pca[:,0])+1) # color map to the number of frames in trajectory 
	plt.scatter(pca[:,0], pca[:,1], marker='x', c=col)
	plt.xlabel('PC1')
	plt.ylabel('PC2')
	plt.title(title)
	cbar = plt.colorbar()
	cbar.set_label('Time [ps]')
	fig.savefig(fname)
	
	## # plot 1 and 3 PC
	fname = out_dir+'/'+file_name+'1_3'+'.png'
	fig=plt.figure()
	plt.scatter(pca[:,0], pca[:,2], marker='x', c=col)
	plt.xlabel('PC1')
	plt.ylabel('PC3')
	plt.title(title)
	cbar = plt.colorbar()
	cbar.set_label('Time [ps]')
	fig.savefig(fname)
	
	## # plot 2 and 3 PC
	fname = out_dir+'/'+file_name+'2_3'+'.png'
	fig=plt.figure()
	plt.scatter(pca[:,1], pca[:,2], marker='x', c=col)
	plt.xlabel('PC2')
	plt.ylabel('PC3')
	plt.title(title)
	cbar = plt.colorbar()
	cbar.set_label('Time [ps]')
	fig.savefig(fname)
	return;

## write PCs 
def write_pcs(file_name, pca, out_dir):
	'write PCs and explained_variance_ratio_. takes name of the file to write and pca object name'
	fname = ''
	fname = out_dir+'/'+file_name+'.agr'
	#print type(pca)
	e_ratio = pca.explained_variance_ratio_
	e_ratio = e_ratio*100   # to make it percent
	
	np.savetxt(fname, e_ratio)
	ef = open(fname, 'r')
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
	
	ef = open(fname, 'w')
	ef.write('#'+title+'\ton\t'+my_time+'\n'+legends+'\n'+ef_cont)
	ef.close()
	return;
	
