#!/usr/bin/env python

import numpy as np
import os
from datetime import datetime
import h5py
import sys
from mpi4py import MPI

from yaff import *

from mylammps import *

from yaff import ForceField, System, NeighborList, Scalings, PairPotLJ,     ForcePartPair, ForcePart, log
log.set_level(log.low)


# Setup MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

# Set random seed, important to get the same velocities for all processes
np.random.seed(5)
# Turn off logging for all processes, it can be turned on for one selected process later on

if __name__=='__main__':
    f_res = h5py.File('restart_5.h5', mode='r')
    system = System.from_hdf5(f_res)
    ff = ForceField.generate(system, 'pars.txt', rcut=15*angstrom, alpha_scale=3.2, gcut_scale=1.5, smooth_ei=True)

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
    pair_lammps = ForcePartLammps(system, pppm_accuracy=1e-6,fn_table='lammps_smoothei2.table', comm=comm, scalings=[0.0,0.0,1.0,0.0,0.0,1.0])

    ff = ForceField(system,[part_valence,pair_lammps,pair_gauss],nlist)

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

    f = h5py.File('traj_6.h5', mode='w')
    hdf = HDF5Writer(f, step=1000)
    g = h5py.File('restart_6.h5',mode='w')
    hdf_restart = RestartWriter(g, step=100000)

    vsl = VerletScreenLog(step=5000)
    md = VerletIntegrator(ff, hooks=[hdf, vsl, hdf_restart], restart_h5 = f_res)
    md.run(2000000)
