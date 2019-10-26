#!/usr/bin/env python

import numpy as np
import os
import h5py
from mpi4py import MPI
from glob import glob
from yaff import *
from mylammps import *
from yaff import ForceField, System, NeighborList, Scalings, PairPotLJ,     ForcePartPair, ForcePart, log
from tailcorr import ForcePartTailCorrection

# Setup MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

# Set random seed, important to get the same velocities for all processes
np.random.seed(5)
# Turn off logging for all processes, it can be turned on for one selected process later on
log.set_level(log.silent)
if rank==0: log.set_level(log.medium)

f_res = h5py.File('restart_1.h5', mode='r')
system = System.from_hdf5(f_res)


fns = []
for fn in os.listdir(os.getcwd()):
    if fn.startswith('pars') and fn.endswith('.txt'):
        fns.append(fn)

rcut = 12*angstrom

ff = ForceField.generate(system, fns, rcut=rcut, smooth_ei=False, gcut_scale=1.5, alpha_scale=3.2)

tailcorr = True
if tailcorr:
    for part in ff.parts:
        if part.name=='pair_mm3':
            sigmas = part.pair_pot.sigmas.copy()
            epsilons = part.pair_pot.epsilons.copy()
            onlypaulis = part.pair_pot.onlypaulis.copy()
            tr = part.pair_pot.get_truncation()
            pair_pot = PairPotMM3(sigmas, epsilons, onlypaulis, rcut, tr=tr)
            part = ForcePartPair(ff.system, ff.nlist, part.scalings, pair_pot)
            part_tail = ForcePartTailCorrection(ff.system, part.pair_pot)
            ff.add_part(part_tail)

# Replace non-covalent part with LAMMPS
for part in ff.parts:
    if part.name=='valence': part_valence = part
    elif part.name=='pair_ei': part_ei = part
# Write LAMMPS system file
write_lammps_data(system)
nlist = BondedNeighborList(system)
tr = part_ei.pair_pot.get_truncation()
if tr is not None:
    pair_pot = PairPotEI(system.charges, 0.0, part_ei.pair_pot.rcut, tr=Switch3(tr.width), radii=system.radii)
else:
    pair_pot = PairPotEI(system.charges, 0.0, part_ei.pair_pot.rcut, tr=None, radii=system.radii)

scalings = Scalings(system, scale1=1.0, scale2=1.0, scale3=0.0)
pair_gauss = ForcePartPair(system,nlist,scalings,pair_pot)
pair_lammps = ForcePartLammps(system, pppm_accuracy=1e-6,fn_table='lammps_smoothei2.table', comm=comm, scalings=[0.0,0.0,1.0,0.0,0.0,1.0],fn_log='none')

ff = ForceField(system,[part_valence,pair_lammps,pair_gauss],nlist)

tailcorr = True
if tailcorr:
    for part in ff.parts:
        if part.name=='pair_mm3':
            sigmas = part.pair_pot.sigmas.copy()
            epsilons = part.pair_pot.epsilons.copy()
            onlypaulis = part.pair_pot.onlypaulis.copy()
            tr = part.pair_pot.get_truncation()
            pair_pot = PairPotMM3(sigmas, epsilons, onlypaulis, rcut, tr=tr)
            part = ForcePartPair(ff.system, ff.nlist, part.scalings, pair_pot)
            part_tail = ForcePartTailCorrection(ff.system, part.pair_pot)
            ff.add_part(part_tail)

# Only write output from process 0
if rank==0:
    f = h5py.File('traj_2.h5', mode='w')
    hdf = HDF5Writer(f, step=1000)
    g = h5py.File('restart_2.h5',mode='w')
    hdf_restart = RestartWriter(g, step=100000)

vsl = VerletScreenLog(step=5000)

# Only write output from process 0
if rank==0:
    log.set_level(log.medium)
    hooks = [hdf,hdf_restart,vsl]
else:
    log.set_level(log.silent)
    hooks = [vsl]

md = VerletIntegrator(ff, hooks=hooks, restart_h5 = f_res)
md.run(500000)
