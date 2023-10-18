#!/usr/bin/env python3
'''
Combines several trajectories before performing PCA
Some PCA material is borrowed from 'http://mdtraj.org/latest/examples/pca.html'
'''
import sys
import argparse
import numpy as np
import mdtraj as md
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
import matplotlib
import seaborn as sns
matplotlib.use('Agg')
from matplotlib import pyplot as plt

__author__ = "Olivier Sheik Amamuddy"
__copyright__ = "Copyright 2019, Research Unit in Bioinformatics"
__license__ = "GNU GPL 3.0"
__version__ = "1.2"
__maintainer__ = "Olivier Sheik Amamuddy"
__email__ = "oliserand@gmail.com"
__date__ = "24th May 2021"

def parse_args():
    """
    The main argument parser
    """
    parser = argparse.ArgumentParser(description="Combines several trajectories\
            and aligns them before performing a single PCA decomposition on the\
            entire array. One important requirement is that the data used\
            should be comparable.")
    parser.add_argument('--topologies', '-s', type=str, nargs='+', required=True,
                        help="A list of topology file names")
    parser.add_argument('--trajectories', '-f', type=str, nargs='+', required=True,
                        help="The corresponding trajectory file names, in the \
                        same order as the topologies")
    parser.add_argument('--selection', type=str,
                        default="all",
                        help="Atom selection (Use MDTraj selection DSL) \
                        (default='all'). Note: By default, only CA \
                        atoms are used. This is usually most applicable when \
                        specifying residue(s) and/or chain ID(s). Selection is \
                        performed after having done a global fitting of the \
                        entire structure.")
    parser.add_argument('--stride', type=int, default=1,
                        help="Step size through each trajectory (default=1)")
    parser.add_argument('--suffix', type=str, default="pcafig",
                        help="File name suffix for the figures")
    parser.add_argument('--n_clusters', type=int, default=3,
                        help="The expected number of protein clusters to\
                        extract.")
    parser.add_argument('--ignn', type=int, default=0,
                        help="The number of N-terminus residues to ignore in \
                        PCA calculations (default=0)")
    parser.add_argument('--ignc', type=int, default=3,
                        help="The number of C-terminus residues to ignore in \
                        PCA calculations (default=3)")
    args = parser.parse_args()
    return args

def check_compatibility(num_frames):
    """Check if trajectories are compatible given this list of frame numbers"""
    assert isinstance(num_frames, list), "A list has to be provided"
    if len(num_frames) > 1:
        if len(set(num_frames)) > 1:
            raise ValueError("Trajectory length mismatch. Ensure all trajectories "
                    "have the same number of frames.")
            sys.exit(1)

def clean_protein_termini(traj, selection="name CA", ignn=0, ignc=0):
    """
    Reduces protein to CA atoms and removes the C and N termini in each chain
     Input
     traj: A trajectory
     selection: The selection to use. The default parameter will usually be the
                best option for getting protein atoms.
     ignn: Number of N-terminus residues to ignore
     ignc: Number of C-terminus residues to ignore
     Returns
     cleaned trajectory
    """
    atoms = traj.top.select(selection)
    traj = traj.atom_slice(atoms)
    chains = traj.top.chains
    atom_indices = []
    for chain in chains:
        atoms = chain.atoms
        index = [a.index for a in atoms]
        index = index[ignn:(len(index)-ignc)]
        atom_indices += index
    traj.atom_slice(atom_indices, inplace=True)
    return traj

def get_kmean_centroids(pcs, traj, outbasename, n_clusters=3,
                        n_init=50, max_iter=1000, marker_size=0.5,
                        marker="o"):
    """
    Uses partitioning, before extracting a defined number of representative\
    conformations for each partition
    """
    coords = pcs[:,:2]
    kmeans = KMeans(n_clusters=n_clusters, n_init=n_init, max_iter=max_iter)
    kmeans.fit(coords)
    centroids = kmeans.cluster_centers_
    # Plot
    plt.scatter(pcs[:, 0], pcs[:, 1], s=marker_size, c=traj.time, cmap="viridis")
    cbar = plt.colorbar()
    cbar.set_label('Time (ps)')
    # Get closest points to each centroid
    for centroid in centroids:
        displacements = coords - centroid
        position = np.linalg.norm(displacements, axis=1).argmin()
        nearest_neighbor = coords[position]
        conformation_time =  traj.time[position]
        plt.text(nearest_neighbor[0], nearest_neighbor[1],
                s=int(conformation_time), c="red", fontsize=8,
                bbox=dict(facecolor='white', pad=0.5, alpha=0.5, lw=0.2))
        plt.scatter(nearest_neighbor[0], nearest_neighbor[1], marker=marker,
                edgecolor="red", facecolors="none")
        # Dump structure
        currtraj = traj[position]
        currtraj.save_pdb("{}_kmeans_centroid_{}ps.pdb".format(outbasename, int(currtraj.time[0])))
    return centroids

def get_energy_basins(pcs, traj, outbasename, marker="o",
                      kde_levels=30, plt=plt):
    """
    Searches the landscape, and plots the lowest energy basins
     Input:
     pcs: Array from PCA
     timelabels: Time (array) obtained from the mdtraj object
     marker: plot character for showing the location of the basin
     kde_levels: Number of levels for 2D landscape
     Returns:
     The basin centroids
    """
    subplot = sns.kdeplot(x=pcs[:, 0], y=pcs[:, 1], shade=True, alpha=0.2,
                          levels=kde_levels, cmap="jet_r")
    basin_paths = subplot.collections[-1].get_paths()
    basins = [basin.vertices for basin in basin_paths]
    basin_centroids = [basin.mean(axis=0) for basin in basins]
    # Search basin for nearest nearest neighbor
    coords = pcs[:,:2]
    for basin_centroid in basin_centroids:
        displacements = coords - basin_centroid
        position = np.linalg.norm(displacements, axis=1).argmin()
        nearest_neighbor = coords[position]
        conformation_time =  traj.time[position]
        plt.text(nearest_neighbor[0], nearest_neighbor[1],
                s=int(conformation_time), c="blue", fontsize=8,
                bbox=dict(facecolor='white', pad=0.5, alpha=0.5, lw=0.2))
        plt.scatter(nearest_neighbor[0], nearest_neighbor[1], marker=marker,
                edgecolor="black", facecolors="none")
        currtraj = traj[position]
        currtraj.save_pdb("{}_basin_{}ps.pdb".format(outbasename, int(currtraj.time[0])))
    return basin_centroids

def plot_graphs(pcs, outbasename, traj, title="Essential dynamics plot",
                xlim=None, ylim=None, percent_variance=None,
                n_clusters=3):
    """
    Generates scatter plots for PCA
    """
    fig = plt.figure(dpi=500)
    plt.xlim(xlim)
    plt.ylim(ylim)
    # Search basins and get KDE samples
    get_energy_basins(pcs, traj=traj, outbasename=outbasename)
    get_kmean_centroids(pcs, traj=traj, outbasename=outbasename,
                    n_clusters=n_clusters)
    plt.title(title)
    plt.xlabel("PC1 ({:3.2f}%)".format(percent_variance[0]))
    plt.ylabel("PC2 ({:3.2f}%)".format(percent_variance[1]))
    plt.tight_layout()
    plt.savefig("{}.png".format(outbasename))

def main(args):
    """
    Program main
    """
    sns.set_style("whitegrid")
    plt.rcParams["font.family"] = 'serif'
    pca = PCA()
    ntrajectories = len(args.trajectories)
    num_frames = []
    all_frames = []
    all_times = []
    topologies = args.topologies
    trajectories = args.trajectories
    n_clusters = args.n_clusters
    # Set a reference frame
    ref_frame = md.load(args.trajectories[0], top=args.topologies[0])
    ref_frame = clean_protein_termini(ref_frame, ignn=args.ignn, ignc=args.ignc)
    # Read and align trajectories
    for i, topology in enumerate(topologies):
        trajectory_filename = trajectories[i]
        traj = md.load(trajectory_filename, top=topology)
        traj = clean_protein_termini(traj, ignn=args.ignn, ignc=args.ignc)
        traj = traj[::args.stride]
        traj = traj.superpose(ref_frame, 0)
        # Post alignment selection
        if args.selection != "all":
            traj = traj.atom_slice(traj.top.select(args.selection))
        xyz_reshaped = traj.xyz.reshape(traj.n_frames, traj.n_atoms*3)
        all_frames.append(xyz_reshaped)
        all_times.append(traj.time)
        num_frames.append(traj.n_frames)
        check_compatibility(num_frames)
    xyz_concat = np.concatenate(all_frames)
    pca.fit(xyz_concat)
    percent_variance = pca.explained_variance_ratio_[:2]*100
    pcs = pca.transform(xyz_concat)
    # Setting plot limits
    xmin = np.min(pcs[:, 0])
    xmax = np.max(pcs[:, 0])
    ymin = np.min(pcs[:, 1])
    ymax = np.max(pcs[:, 1])
    # Iterate in chunks of trajectory length (time)
    increment = int(pcs.shape[0]/ntrajectories)
    for idx, start in enumerate(range(0, pcs.shape[0], increment)):
        end = start+increment
        # Slice PCs for current trajectory
        pcs_currtraj = pcs[start:end, :]
        # Get cluster representatives
        outbasename = args.trajectories[idx][:-4]
        # Reload full trajectory
        traj = md.load(args.trajectories[idx],
                       top=args.topologies[idx])[::args.stride]
        # Plotting
        plot_graphs(pcs=pcs_currtraj, traj=traj, outbasename=outbasename,
                    xlim=(xmin, xmax), ylim=(ymin, ymax),
                    percent_variance=percent_variance, n_clusters=n_clusters)

if __name__ == "__main__":
    ARGS = parse_args()
    main(ARGS)
