#! /usr/bin/env python

from molmod.units import *
from yaff import *
import numpy as np
import h5py

f_res = h5py.File('restart_1.h5', mode='r')
system = System.from_hdf5(f_res)
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

f = h5py.File('traj_2.h5', mode='w')
hdf = HDF5Writer(f, step=100)
g = h5py.File('restart_2.h5', mode='w')
hdf_restart = RestartWriter(g, step=50000)

vsl = VerletScreenLog(step=5000)
md = VerletIntegrator(ff, hooks=[hdf, vsl, hdf_restart], restart_h5 = f_res)
md.run(400000)
