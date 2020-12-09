#! /usr/bin/env python

from molmod.units import *
from yaff import *
import numpy as np
import h5py


f = h5py.File('traj_1.h5', mode='r')
start = 1
end = -1
step = 1
volume = np.array(f['/trajectory/volume'][start:end:step])/(angstrom**3)
g = open('volume.txt','w')
for i in xrange(len(volume)):
    g.write(str(volume[i])+'\n')
