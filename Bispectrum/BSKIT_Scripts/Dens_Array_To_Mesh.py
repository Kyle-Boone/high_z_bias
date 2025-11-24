import os
import argparse
import numpy as np
from nbodykit.lab import ArrayMesh
from abacusnbody.analysis.tsc import tsc_parallel

parser = argparse.ArgumentParser()
parser.add_argument('dens_file',
                    help='File of positional galaxy data')
parser.add_argument('write_file',
                    help='File to write the density grid to temporarily')
parser.add_argument('boxsize', type = float,
                    help='BoxSize in Mpc/h')
args = parser.parse_args()

dens_file = args.dens_file
write_file = args.write_file
Lbox = args.boxsize

d = np.load(dens_file)
d = d - np.average(d)

# Create the MeshSource object
mesh = ArrayMesh(d, BoxSize=Lbox)

mesh.save(write_file)
