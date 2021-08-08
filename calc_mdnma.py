#!/usr/bin/env python3
"""
Computes ANM from MD. Writes the mapping file for NGL and can display 
the mode animation using a multi-PDB file
"""
import argparse
import mdtraj as md
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

__author__ = "Olivier Sheik Amamuddy"
__license__ = "GNU GPL 3.0"
__version__ = "1.0.0"
__maintainer__ = "Olivier Sheik Amamuddy"
__email__ = "oliserand@gmail.com"
__status__ = "Development"
__date__ = "22nd February 2021"

COLORS = {0:np.array([255,0,0])/255, 1:np.array([0,0,255])/255,
          2:np.array([204,204,0])/255, 3:np.array([204,0,204])/255,
          4:np.array([255,255,153])/255, 5:np.array([255,51,51])/255,
          6:np.array([0,255,255])/255, 7:np.array([255,204,229])/255,
          8:np.array([224,224,224])/255, 9:np.array([255,102,255])/255,
          10:np.array([153,153,0])/255, 11:np.array([102,178,255])/255,
          12:np.array([153,255,255])/255, 13:np.array([153,204,255])/255,
          14:np.array([255,255,204])/255, 15:np.array([128,255,0])/255,
          16:np.array([204,255,153])/255, 17:np.array([204,0,204])/255,
          18:np.array([255,51,255])/255, 19:np.array([255,204,255])/255}

def arg_parser():
    """Argument parser"""
    parser = argparse.ArgumentParser(
            description="Calculates normal modes from PCA")
    parser.add_argument("topology", type=str,
                         help="MD topology filename")
    parser.add_argument("trajectory", type=str,
                        help="MD trajectory filename")
    parser.add_argument("--mode", type=int, default=1,
            help="Number of residues to skip at the C-terminus (default:1)")
    parser.add_argument("--selection", type=str, default="all",
            help="Residue selection using MDTraj syntax (default:name CA)")
    parser.add_argument("--ignn", type=int, default=0,
                        help="Number of residues to skip at the N-terminus\
                             (default:0)")
    parser.add_argument("--ignc", type=int, default=3,
            help="Number of residues to skip at the C-terminus (default:3)")
    parser.add_argument("--arrow_scale", type=float, default=6.0,
            help="Factor by which to scale each arrow in NGLView (default:6)")
    parser.add_argument("--outfilename", type=str,
                        default="NGL_ARROWS.js",
                        help="Output filename for NGLView")
    parser.add_argument("--animate", action="store_true",
            help="Generate an animation (PDB format) by adding the modes \
                    (default:no)")
    return parser.parse_args()

def write_header(fobject, topology):
    """Writes header"""
    header = """$.when(
$.get("/api/mdtask/jobs/" + job.selected_job().JobID() + "/topomdnma", function(responseTopo) { }),
).then(function(responseTopo) {
var stage = new NGL.Stage( "viewport_mdnma_struc" );
window.addEventListener( "resize", function( event ){
stage.handleResize();
}, false );
var shape = new NGL.Shape("shape", { disableImpostor: true, radialSegments: 10 });
stage.loadFile( "%s" ).then( function( o ){
o.addRepresentation( "trace" );
o.autoView();
} );""" % topology
    fobject.write(header)

def write_footer(fobject):
    """Writes footer"""
    footer = """var shapeComp = stage.addComponentFromObject(shape)
shapeComp.addRepresentation("buffer", { wireframe: false })
var stringBlobTopo = new Blob( [ responseTopo ], { type: "text/plain"} );
stage.loadFile( stringBlobTopo, { asTrajectory: true, ext: "pdb" } ).then(function (o) {
o.autoView()
})
stage.setParameters( { backgroundColor: "white", hoverTimeout: -1 } );
});\n"""
    fobject.write(footer)

def check_mode(mode, lmax):
    """Checks mode validity"""
    if mode < 0:
        raise IndexError("ERROR: Mode cannot be less than 1")
    if mode >= lmax:
        raise IndexError("ERROR: Mode number exceeded")

def animate_nma(evects, mdtraj_top, mode=0, scalar=0.02, nframes=50, maxrmsd=0.1):
    """Animates NMA by writing a multi-PDB file"""
    evect = evects[:, mode]
    first_frame = mdtraj_top
    trajectory = []
    ev_reshaped = evect.reshape(evect.shape[0]//3, 3)
    atoms = first_frame.top.select("all")
    for i in range(nframes):
        frame = first_frame.atom_slice(atoms)
        frame.xyz += i*scalar*ev_reshaped
        rmsd = md.rmsd(frame, first_frame)
        if rmsd[0] < maxrmsd:
            trajectory.append(frame)
        else:
            break
    frames = first_frame.join(trajectory)
    frames.save_pdb("pcamode_{}.pdb".format(mode+1))
    print("INFO: Wrote file \"pcamode_{}.pdb\"".format(mode+1))

def calc_nma(trajectory, topology, selection="all",
             ignn=0, ignc=0, mode=0, animate=True):
    """Calculate NMA"""
    traj = md.load(trajectory, top=topology)
    atoms = traj.top.select("name CA")
    atoms = atoms[ignn:len(atoms)-ignc]
    traj.atom_slice(atoms, inplace=True)
    traj.atom_slice(traj.top.select(selection), inplace=True)
    check_mode(mode, traj.n_residues)
    traj = traj.superpose(reference=traj)
    atoms = [a.index for a in traj.top.atoms]
    cartesian = traj.xyz.reshape(traj.n_frames, traj.n_residues*3).transpose()
    cartesian_centered = cartesian - \
            cartesian.mean(axis=1).reshape(cartesian.shape[0], 1)
    covmat = cartesian_centered.dot(cartesian_centered.transpose())
    evals_, evects_ = np.linalg.eigh(covmat)
    order = np.argsort(evals_)[::-1]
    evals = evals_[order]
    evects = evects_[:,order]
    perc_exp_var_ = 100*(evals/evals.sum())
    perc_exp_var = perc_exp_var_[mode].round(2)
    evect = evects[:, mode]
    _, axis = plt.subplots(figsize=(8, 3), dpi=450)
    axis.bar(x=range(50), height=perc_exp_var_[:50])
    axis.set_xlabel("Indices for the first 50 modes")
    axis.set_ylabel("Percentage explained variance")
    plt.tight_layout()
    plt.savefig("pca_eigenvalues.png".format(mode+1))
    print("INFO: Wrote file \"pca_eigenvalues.png\"".format(mode+1))
    if animate:
        animate_nma(evects, traj[0], mode=mode)
    return evect

def write_ngl(nma, args):
    """Write ngl file for the mapping"""
    protein = md.load(args.topology)
    protein = protein.atom_slice(protein.top.select("name CA"))
    protein = protein.atom_slice(
            protein.top.select("resid >= {0} and resid < {1} ".format(
                args.ignn, protein.n_residues-args.ignc)))
    topology = protein.atom_slice(protein.top.select(args.selection))
    arrow_starts = topology.xyz[0,:,:]*10
    arrow_vects = nma.reshape(nma.shape[0]//3, 3)*10
    arrow_ends = arrow_starts + arrow_vects*args.arrow_scale
    data, _ = topology.top.to_dataframe()
    chains = data["chainID"].values
    with open(args.outfilename, "w") as _:
        write_header(_, args.topology)
        text = "\n"
        for i in range(arrow_starts.shape[0]):
            text += "shape.addArrow([{0},{1},{2}], [{3},{4},{5}], [{6},{7},{8}], {9})\n".format(
                    float(arrow_starts[i,0]), float(arrow_starts[i,1]), float(arrow_starts[i,2]),
                    float(arrow_ends[i,0]), float(arrow_ends[i,1]), float(arrow_ends[i,2]),
                    COLORS[chains[i]][0],COLORS[chains[i]][1],COLORS[chains[i]][2], 0.5)
        _.write(text)
        write_footer(_)
        print("INFO: Wrote file \"{}\"".format(args.outfilename))

def main(args):
    """Main program"""
    nma = calc_nma(args.trajectory, args.topology,
                   selection=args.selection, ignn=args.ignn,
                   ignc=args.ignc, mode=args.mode-1,
                   animate=args.animate)
    write_ngl(nma, args)

if __name__ == "__main__":
    plt.rcParams["font.family"] = "serif"
    args = arg_parser()
    main(args)
