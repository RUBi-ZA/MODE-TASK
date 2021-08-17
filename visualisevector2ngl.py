#!/usr/bin/env python3
""" Library/ script for converting TCL arrows to NGL. Mainly useful for embedded NGLview code
Olivier Sheik Amamuddy
18th Feb 2021
"""
import re 
import argparse
import numpy as np
from subprocess import call

colors = {"red":np.array([255,0,0])/255, "blue":np.array([0,0,255])/255,
          "ochre":np.array([204,204,0])/255, "purple":np.array([204,0,204])/255,
           "yellow":np.array([255,255,153])/255, "red":np.array([255,51,51])/255,
           "cyan":np.array([0,255,255])/255, "pink":np.array([255,204,229])/255,
           "silver":np.array([224,224,224])/255, "violet":np.array([255,102,255])/255,
           "ochre":np.array([153,153,0])/255, "blue2":np.array([102,178,255])/255,
           "cyan2":np.array([153,255,255])/255, "iceblue":np.array([153,204,255])/255,
           "lime":np.array([255,255,204])/255, "green2":np.array([128,255,0])/255,
           "green3":np.array([204,255,153])/255, "violet":np.array([204,0,204])/255,
           "violet2":np.array([255,51,255])/255, "mauve":np.array([255,204,255])/255} 


def get_arrow_string(text, color=None):
    text = text.replace("draw arrow ", "shape.addArrow(")
    text = text.replace(" ",",").replace("{","[").replace("}","]").replace("\n","")
    text = "{},[{},{},{}], 1)\n".format(text, *color)
    return text

def write_header(fobject):
    header = """$.when(
$.get("/api/mdtask/jobs/" + job.selected_job().JobID() + "/toponmaa", function(responseTopo) { }),
).then(function(responseTopo) {
var stage = new NGL.Stage( "viewport_vm_out" );
window.addEventListener( "resize", function( event ){
stage.handleResize();
}, false );
var shape = new NGL.Shape("shape", { disableImpostor: true, radialSegments: 10 })
"""
    fobject.write(header)

def write_footer(fobject):
    footer = """var shapeComp = stage.addComponentFromObject(shape)
shapeComp.addRepresentation("buffer", { wireframe: false })
var stringBlobTopo = new Blob( [ responseTopo ], { type: "text/plain"} );
stage.loadFile( stringBlobTopo, { asTrajectory: true, ext: "pdb" } ).then(function (o) {
o.autoView()
})
stage.setParameters( { backgroundColor: "white", hoverTimeout: -1 } );
});\n"""
    fobject.write(footer)

def main(args):
    with open(args.filename, "r") as f:
        lines = f.readlines()
        outfile = open(args.outfilename, "w")
        write_header(outfile)
        for line in lines:
            if line.startswith("draw"):
                tokens = line.split()
                if tokens[1] == "color":
                    color = line.split()[-1]
                    color = colors[color]
                elif tokens[1] == "arrow":
                    arrow_text = get_arrow_string(line, color=color)
                    outfile.write(arrow_text)
        write_footer(outfile)
        print("INFO: Wrote file \"{}\"".format(args.outfilename))

def arg_parser():
    parser = argparse.ArgumentParser("Prepares NGL js file from MODE task ")
    parser.add_argument("filename", type=str,
                        help="File from MODE-TASK's visualiseVector.py")
    parser.add_argument("--outfilename", type=str,
                        default="NGL_ARROWS.js",
                        help="Output filename for NGLView")

    args = parser.parse_args()
    return args

if __name__ == "__main__":
    ARGS = arg_parser()
    main(ARGS)
