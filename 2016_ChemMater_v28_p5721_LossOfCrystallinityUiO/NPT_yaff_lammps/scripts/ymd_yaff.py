#! /usr/bin/env python

from molmod.units import *
from yaff import *
import numpy as np
import h5py

system = System.from_file('init.chk')
ff = ForceField.generate(system, 'pars.txt', rcut=15*angstrom, alpha_scale=3.2, gcut_scale=1.5, smooth_ei=True)

ff.part_valence.dlist.forward()
ff.part_valence.iclist.forward()

for I in xrange(ff.part_valence.iclist.nic):
    if ff.part_valence.iclist.ictab[I]['kind'] == 10:
        if ff.part_valence.iclist.ictab[I]['value'] < 0.0:        
            ff.part_valence.iclist.ictab[I]['sign0'] = ff.part_valence.iclist.ictab[I]['sign0']*(-1)
            print 'Found a negative one, adjust!'
            
ff.part_valence.dlist.forward()
ff.part_valence.iclist.forward()
ff.part_valence.vlist.forward()

dof = StrainCellDOF(ff, gpos_rms=1e-8, dpos_rms=1e-6, gcell_rms=1e-8, dcell_rms=1e-6, do_frozen=False)

f = h5py.File('traj_1.h5', mode='w')
hdf = HDF5Writer(f, step=100)
g = h5py.File('restart_1.h5',mode='w')
hdf_restart = RestartWriter(g, step=100000)

temp = 300*kelvin
timestep = 0.5*femtosecond
press=1e5*pascal

thermo = NHCThermostat(temp, timecon=100.0*femtosecond)
baro = MTKBarostat(ff, temp, press, timecon=1000.0*femtosecond)
TBC = TBCombination(thermo, baro)

vsl = VerletScreenLog(step=5000)
md = VerletIntegrator(ff, timestep, hooks=[hdf, TBC, vsl, hdf_restart])
md.run(700000)
