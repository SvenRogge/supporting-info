#!/usr/bin/env python

import numpy as np
import os
import h5py
from yaff import *


log.set_level(log.medium)

system = System.from_file('init.chk')
fns = []
for fn in os.listdir(os.getcwd()):
    if fn.startswith('pars') and fn.endswith('.txt'):
        fns.append(fn)

ff = ForceField.generate(system, fns, rcut=15.0*angstrom, alpha_scale=3.2, gcut_scale=1.5, smooth_ei=True)

f = h5py.File('traj_1.h5', mode='w')
hdf = HDF5Writer(f, step=1000)
g = h5py.File('restart_1.h5',mode='w')
hdf_restart = RestartWriter(g, step=100000)

temp = 300*kelvin
timestep = 0.5*femtosecond
press=1e6*pascal

thermo = NHCThermostat(temp, timecon=100.0*femtosecond)
baro = MTKBarostat(ff, temp, press, timecon=1000.0*femtosecond)
TBC = TBCombination(thermo, baro)
vsl = VerletScreenLog(step=5000)

hooks = [hdf,hdf_restart,TBC,vsl]

md = VerletIntegrator(ff, timestep, hooks=hooks)
md.run(500000)
