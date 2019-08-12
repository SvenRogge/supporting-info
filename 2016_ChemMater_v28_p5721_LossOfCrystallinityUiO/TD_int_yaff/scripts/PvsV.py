#! /usr/bin/env python

from molmod.units import *
from yaff import *
import numpy as np
import h5py

eq_time = 50 # in ps
timestep = 0.05 # in ps

f = h5py.File('traj_conc12.h5', mode='r')
start, end, step = get_slice(f)
start = int(eq_time / timestep)
volume = np.array(f['/trajectory/volume'][-1])/(angstrom**3)
press = np.mean(np.array(f['/trajectory/press'][start:end:step]))/(1e6*pascal)
g = open('../PvsV_eq50ps.csv','a')
g.write(str(volume) + '\t' + str(press) + '\n')
g.close()
