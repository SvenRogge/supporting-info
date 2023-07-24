#!/usr/bin/env python

import numpy as np
import os
import h5py
from mpi4py import MPI
from glob import glob
from yaff import *
from molmod import kjmol

# Setup MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

# Set random seed, important to get the same velocities for all processes
np.random.seed(5)
# Turn off logging for all processes, it can be turned on for one selected process later on
log.set_level(log.silent)
if rank==0: log.set_level(log.medium)

system = System.from_file('ZIF8_6SF6.chk')

rcut = 12*angstrom
fns = 'pars.txt'

ff_yaff = ForceField.generate(system, fns, rcut=rcut, smooth_ei=False, gcut_scale=1.5, alpha_scale=3.2, tailcorrections=True)

ff = swap_noncovalent_lammps(ff_yaff, fn_system='lammps.dat', fn_log="log.lammps", suffix='', fn_table='lammps_smoothei2.table', comm=comm)

print('Yaff energy')
print(ff_yaff.compute()/kjmol)

print('LAMMPS energy')
print(ff.compute()/kjmol)

