#!/usr/bin/env python

import numpy as np
import os
import h5py
from mpi4py import MPI
from glob import glob

from yaff import *

from mylammps import *

from yaff import ForceField, System, NeighborList, Scalings, PairPotLJ,     ForcePartPair, ForcePart, log
log.set_level(log.low)

from molmod.units import *

# Setup MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

# Set random seed, important to get the same velocities for all processes
np.random.seed(5)
# Turn off logging for all processes, it can be turned on for one selected process later on
log.set_level(log.silent)

if __name__=='__main__':
    system = System.from_file('init.chk')
    log.set_level(log.silent)
    ff = ForceField.generate(system, 'pars.txt', rcut=15*angstrom, alpha_scale=3.2, gcut_scale=1.5, smooth_ei=True)

    print ff.compute()/kjmol

    # Replace non-covalent part with LAMMPS
    for part in ff.parts:
        if part.name=='valence': part_valence = part
        elif part.name=='pair_ei': part_ei = part
    # Write LAMMPS system file
    write_lammps_data(system)
    nlist = BondedNeighborList(system)
    pair_pot = PairPotEI(system.charges, 0.0, part_ei.pair_pot.rcut, tr=Switch3(part_ei.pair_pot.get_truncation().width), radii=system.radii)
    scalings = Scalings(system, scale1=1.0, scale2=1.0, scale3=0.0)
    pair_gauss = ForcePartPair(system,nlist,scalings,pair_pot)
    pair_lammps = ForcePartLammps(system, pppm_accuracy=1e-6,fn_table='lammps_smoothei2.table', comm=comm, scalings=[0.0,0.0,1.0,0.0,0.0,1.0],fn_log='none')

    ff = ForceField(system,[part_valence,pair_lammps,pair_gauss],nlist)

    print ff.compute()/kjmol
    for part in ff.parts:
        print "%20s %20.12f" % (part.name,part.energy/kjmol)

