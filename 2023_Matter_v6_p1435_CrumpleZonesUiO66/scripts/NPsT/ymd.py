#!/usr/bin/env python

import numpy as np
import os
import h5py
from mpi4py import MPI
from glob import glob
from yaff import *

# Setup MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

# Set random seed, important to get the same velocities for all processes
np.random.seed(5)
# Turn off logging for all processes, it can be turned on for one selected process later on
log.set_level(log.silent)
if rank==0: log.set_level(log.medium)

sys = System.from_file('init.chk')
system = sys.supercell(2,2,2)

fns = []
for fn in os.listdir(os.getcwd()):
    if fn.startswith('pars') and fn.endswith('.txt'):
        fns.append(fn)

rcut = 12*angstrom

ff_yaff = ForceField.generate(system, fns, rcut=rcut, smooth_ei=False, gcut_scale=1.5, alpha_scale=3.2, tailcorrections=True)
ff = swap_noncovalent_lammps(ff_yaff, fn_system='lammps.dat', fn_log="log.lammps", suffix='', fn_table='lammps_smoothei2.table', comm=comm)

# Only write output from process 0
if rank==0:
    f = h5py.File('traj_1.h5', mode='w')
    hdf = HDF5Writer(f, step=1000)
    g = h5py.File('restart_1.h5',mode='w')
    hdf_restart = RestartWriter(g, step=100000)

temp = 300*kelvin
timestep = 0.5*femtosecond
press=0e6*pascal

thermo = NHCThermostat(temp, timecon=100.0*femtosecond)
baro = MTKBarostat(ff, temp, press, timecon=1000.0*femtosecond, vol_constraint=False)
TBC = TBCombination(thermo, baro)
vsl = VerletScreenLog(step=5000)

# Only write output from process 0
if rank==0:
    log.set_level(log.medium)
    hooks = [hdf,hdf_restart,TBC,vsl]
else:
    log.set_level(log.silent)
    hooks = [TBC,vsl]

md = VerletIntegrator(ff, timestep, hooks=hooks)
md.run(1000000)
