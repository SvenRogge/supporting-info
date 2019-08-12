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


# Setup MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

# Set random seed, important to get the same velocities for all processes
np.random.seed(5)
# Turn off logging for all processes, it can be turned on for one selected process later on
log.set_level(log.silent)
if rank==0: log.set_level(log.medium)

if __name__=='__main__':
    system = System.from_file('init.chk')
    fns = []
    for fn in os.listdir(os.getcwd()):
        if fn.startswith('pars') and fn.endswith('.txt'):
            fns.append(fn)
    ff = ForceField.generate(system, fns, rcut=15.0*angstrom, alpha_scale=3.2, gcut_scale=1.5, smooth_ei=True)

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

    dof = StrainCellDOF(ff, gpos_rms=1e-8, dpos_rms=1e-6, gcell_rms=1e-8, dcell_rms=1e-6, do_frozen=False)

    # Only write output from process 0
    if rank==0:
        f = h5py.File('traj_1.h5', mode='w')
        hdf = HDF5Writer(f, step=1000, flush=10)
        g = h5py.File('restart_1.h5',mode='w')
        hdf_restart = RestartWriter(g, step=100000,flush=1)

    temp = 300*kelvin
    timestep = 0.5*femtosecond
    press=1e6*pascal

    thermo = NHCThermostat(temp, timecon=100.0*femtosecond)
    baro = MTKBarostat(ff, temp, press, timecon=1000.0*femtosecond, vol_constraint=True)
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
    md.run(1200000)
