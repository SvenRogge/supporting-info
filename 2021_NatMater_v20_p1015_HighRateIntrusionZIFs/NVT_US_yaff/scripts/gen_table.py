#!/usr/bin/env python

import os, sys, h5py
import numpy as np
import h5py as h5
from mpi4py import MPI
from glob import glob

#from yaff import System, ForceField, HDF5Writer, RestartWriter, XYZWriter,\
#    NHCThermostat, MTKBarostat, TBCombination, VerletScreenLog, VerletIntegrator, log,\
#    ForcePart, PairPotLJ, PairPotMM3, ForcePartPair, swap_noncovalent_lammps, write_lammps_table
from yaff import *
from molmod.units import *
from molmod.constants import *
from ghostatoms import GhostForceField

def get_ff(rcut=12.0*angstrom, dotail=True):
    # Load initial system
    system = System.from_file('init.chk')
    # Construct force field
    ff_fns = ['pars.txt']
    ff = ForceField.generate(system, ff_fns, rcut=rcut, alpha_scale=3.2, gcut_scale=1.5, tailcorrections=dotail)
    # Get the indexes of the ghost atoms (M atoms of TIP4P water)
    ghost_indexes = np.array([iatom for iatom in range(system.natom) if system.get_ffatype(iatom)=='TM'])
    ghostff = GhostForceField(ff, ghost_indexes, write_ghost_pos, write_ghost_gpos)
    return ghostff

def write_ghost_pos(system, ghost_indexes):
    '''
    Update M site position of TIP4P water based on position of other atoms in molecule
    Assumes that atoms are always ordered as O-H-H-M

        r_M = r_O + d_OM^rel/2 * [ (1+d02/d01)*r01 + (1+d01/d02)*r02 ]
    '''
    d_om_rel = 0.13194
    for iatom in ghost_indexes:
        # Vector pointing from O to H1
        r01 = system.pos[iatom-1] - system.pos[iatom-3]
        system.cell.mic(r01)
        d01 = np.linalg.norm(r01)
        # Vector pointing from O to H2
        r02 = system.pos[iatom-2] - system.pos[iatom-3]
        system.cell.mic(r02)
        d02 = np.linalg.norm(r02)
        # Set M position
        system.pos[iatom] = system.pos[iatom-3] + 0.5*d_om_rel*((1.0+d02/d01)*r01 + (1.0+d01/d02)*r02)

def write_ghost_gpos(gpos, vtens, system, ghost_indexes):
    d_om_rel = 0.13194
    for iatom in ghost_indexes:
        # Vector pointing from O to H1
        r01 = system.pos[iatom-1] - system.pos[iatom-3]
        system.cell.mic(r01)
        d01 = np.linalg.norm(r01)
        # Vector pointing from O to H2
        r02 = system.pos[iatom-2] - system.pos[iatom-3]
        system.cell.mic(r02)
        d02 = np.linalg.norm(r02)
        # Partial derivatives of M positions
        pdiff_01 = gpos[iatom,:]*(1.0+d02/d01) - r01*np.dot(gpos[iatom,:],(d02/d01/d01/d01*r01-r02/d01/d02))
        pdiff_02 = gpos[iatom,:]*(1.0+d01/d02) - r02*np.dot(gpos[iatom,:],(d01/d02/d02/d02*r02-r01/d02/d01))
        # Apply chain rule
        gpos[iatom-3,:] += gpos[iatom,:]
        gpos[iatom-3,:] -= 0.5*d_om_rel*pdiff_01
        gpos[iatom-3,:] -= 0.5*d_om_rel*pdiff_02
        gpos[iatom-1,:] += 0.5*d_om_rel*pdiff_01
        gpos[iatom-2,:] += 0.5*d_om_rel*pdiff_02
        if vtens is not None:
            r_mo = 0.5*d_om_rel*((1.0+d02/d01)*r01 + (1.0+d01/d02)*r02)
            vtens[:] -= np.outer(gpos[iatom],r_mo)
            vtens[:] += np.outer(0.5*d_om_rel*pdiff_01,r01)
            vtens[:] += np.outer(0.5*d_om_rel*pdiff_02,r02)


if __name__=='__main__':

    # Setup MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    # Set random seed, important to get the same velocities for all processes
    np.random.seed(1)

    # Turn off logging for all processes, it can be turned on for one selected process later on
    log.set_level(log.silent)
    if rank==0: log.set_level(log.medium)

    # General settings
    rcut = 12.0*angstrom
    dotail = True

    # Load initial system
    system = System.from_file('init.chk')
    # Construct force field
    ff = ForceField.generate(system, 'pars.txt', rcut=rcut, alpha_scale=3.2, gcut_scale=1.5, tailcorrections=dotail)
    ff = swap_noncovalent_lammps(ff, fn_system='lammps.dat', fn_log="log.lammps", suffix='', fn_table='rcut_12.0.table', overwrite_table=False, comm=comm)
    # Get the indexes of the ghost atoms (M atoms of TIP4P water)
    ghost_indexes = np.array([iatom for iatom in range(system.natom) if system.get_ffatype(iatom)=='TM'])

    #ff_lammps = swap_noncovalent_lammps(ff, fn_system='lammps.dat', fn_log="log.lammps", suffix='', fn_table='rcut_12.0.table', overwrite_table=False, comm=comm)
    ghostff = GhostForceField(ff, ghost_indexes, write_ghost_pos, write_ghost_gpos) 
