#! /usr/bin/env python

from molmod.units import *
from yaff import *
import numpy as np
import h5py

eq_time = 100 # in ps
timestep = 0.5 # in ps

f1 = h5py.File('traj_1.h5', mode='r')
start, end, step = get_slice(f1)
start = int(eq_time / timestep)

volume1 = np.array(f1['/trajectory/volume'][start:-1:step])/(angstrom**3)
press1 = np.array(f1['/trajectory/press'][start:-1:step])/(1e6*pascal)

volume = np.mean(volume1)
press = np.mean(press1)

g = open('../PvsV_eq100ps_av900ps.csv','a')
g.write(str(volume) + '\t' + str(press) + '\n')
g.close()
