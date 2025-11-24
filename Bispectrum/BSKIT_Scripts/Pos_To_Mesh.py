import os
import h5py
import argparse
import numpy as np
from nbodykit.lab import ArrayMesh
from abacusnbody.analysis.tsc import tsc_parallel

parser = argparse.ArgumentParser()
parser.add_argument('pos_file',
                    help='File of positional galaxy data')
parser.add_argument('write_file',
                    help='File to write the density grid to temporarily')
parser.add_argument('boxsize', type = float,
                    help='BoxSize in Mpc/h')
parser.add_argument('ngrid', type = int,
                    help='N that TSC was done for')
args = parser.parse_args()

pos_file = args.pos_file
write_file = args.write_file
Lbox = args.boxsize
Nmesh = args.ngrid
    
with h5py.File(pos_file, 'r') as f:
    dset = f['dataset']
    pos = dset[:].T
    
d = tsc_parallel(pos.astype(np.float32), Nmesh, Lbox, nthread=32)
d = d.astype(np.float32)
ave = np.average(d)
d = (d / np.average(d)) - 1
d = d - np.average(d)
    
# Create the MeshSource object
mesh = ArrayMesh(d, BoxSize=Lbox)

idx = write_file.rfind('/')+1
data_dir = write_file[:idx]

if not os.path.exists(data_dir):
    os.mkdir(data_dir)
    
mesh.save(write_file)
