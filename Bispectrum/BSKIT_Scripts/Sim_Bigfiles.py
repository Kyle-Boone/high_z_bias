import os
import numpy as np
from nbodykit.lab import ArrayMesh


L = 2000
matter_dir = '/mnt/nvme1/kboone/Sims/Matter_Data/'
bigfile_dir = '/mnt/nvme1/kboone/Sims/BigFiles/'

matter_files = os.listdir(matter_dir)

for matter_file in matter_files:
    ind = matter_file.find('.')
    write_file = bigfile_dir + matter_file[:ind] + '.bigfile'

    d = np.load(matter_dir + matter_file)
    d = d - np.average(d)

    mesh = ArrayMesh(d, BoxSize=L)
    mesh.save(write_file)