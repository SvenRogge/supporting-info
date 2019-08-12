#! /usr/bin/env python

from molmod.units import *
from yaff import *
import numpy as np
import h5py


sys0 = System.from_file('init.chk')
system = sys0.supercell(1,2,1)
ff = ForceField.generate(system, 'pars.txt', rcut=15*angstrom, alpha_scale=3.2, gcut_scale=1.5, smooth_ei=True)
dof = StrainCellDOF(ff, gpos_rms=1e-8, dpos_rms=1e-6, gcell_rms=1e-8, dcell_rms=1e-6, do_frozen=False)

f = h5py.File('traj.h5', mode='w') 
hdf = HDF5Writer(f, step=100)

temp = 300*kelvin
timestep = 0.5*femtosecond
press = 1e5*pascal

thermo = BerendsenThermostat(temp, timecon=100.0*femtosecond)
baro = BerendsenBarostat(ff, temp, press, timecon=10000.0*femtosecond)
TBC = TBCombination(thermo, baro)

vsl = VerletScreenLog(step=5000)
md = VerletIntegrator(ff, timestep, hooks=[hdf, TBC, vsl])
md.run(1600000)
