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
from liblammps import *

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


class ForcePartUpperSphere(ForcePart):
    '''
       ForcePart of a constraining sphere.
    '''
    # TODO: Implement PBCs
    def __init__(self, system, indices, K, r0, rv):
        '''
           **Arguments:**
           system
                An instance of the ``System`` class.
           indices
                Numpy array of atom indices.
           K
                Force constaint of the bias potential.
           r0
                Center of the constraining sphere.
           rv
                Radius of the constraining sphere.
        '''
        ForcePart.__init__(self, 'upperwall_{:.1f}_{:.1f}_{:.1f}'.format(*r0), system)
        self.system = system
        self.indices = indices
        self.K = K
        self.r0 = r0
        self.rv = rv
        assert(len(r0)==3)

        if log.do_medium:
            with log.section('FPINIT'):
                log('Force part: %s' % self.name)
                log.hline()

    def _internal_compute(self, gpos, vtens):
        x = self.system.pos[self.indices, :] - self.r0
        norm = np.linalg.norm(x, axis=1)
        e = 0.5*self.K*np.sum((norm - self.rv)**2 * (norm > self.rv))

        if not (norm > 1e-10).all():
            norm[np.where(norm <= 1e-10)[0]] = 1.     

        x *= np.expand_dims(norm > self.rv, axis=1)
        my_gpos =  self.K * (1 - self.rv/np.expand_dims(norm, axis=1)) * x       
                 
        if gpos is not None:
            gpos[self.indices, :] += my_gpos
        if vtens is not None:
            vtens[:] += np.sum(np.expand_dims(x, axis=2)*np.expand_dims(my_gpos, axis=1), axis=0) 
        
        return e


class CVOopDist(CollectiveVariable):
    '''Out-of-plane distance.'''
    def __init__(self, system, indices, indices_ref, indices_plane):
        '''
           **Arguments:**
           system
                An instance of the ``System`` class.
        '''
        self.indices = indices
        self.indices_ref = indices_ref
        self.indices_plane = indices_plane
        CollectiveVariable.__init__(self, 'CVOopDist', system)

    def get_conversion(self):
        return log.length.conversion

    def compute(self, gpos=None, vtens=None):
        # Normal of the plane
        r1 = self.system.pos[self.indices_plane[1]] - self.system.pos[self.indices_plane[0]]
        r2 = self.system.pos[self.indices_plane[2]] - self.system.pos[self.indices_plane[0]]
        #self.system.cell.mic(r1)
        #self.system.cell.mic(r2)
        n = np.cross(r1, r2)
        n /= np.linalg.norm(n)

        # CV
        r = np.mean(self.system.pos[self.indices] , axis=0) - np.mean(self.system.pos[self.indices_ref], axis=0)
        #self.system.cell.mic(r)
        self.value = np.dot(r, n)

        if gpos is not None:
            gpos[:] = 0.0
            gpos[self.indices] = n/len(self.indices)

        if vtens is not None:
            vtens[:] = np.outer(r, n) + np.outer(n, r)

        return self.value


def set_init_config_US(ff, CV0, indices_H2O):

    r1 = ff.system.pos[125] - ff.system.pos[121]
    r2 = ff.system.pos[129] - ff.system.pos[121]
    n = np.cross(r1, r2)
    n /= np.linalg.norm(n)

    pos_H2O = np.mean(ff.system.pos[indices_H2O], axis=0)
    pos_6MR = np.mean(ff.system.pos[np.array([121, 125, 129, 135, 139, 143])], axis=0)

    pos = ff.system.pos.copy()
    pos[indices_H2O] += -pos_H2O + pos_6MR + CV0*n
    ff.update_pos(pos)


class RemoveComMoment(VerletHook):
    def __init__(self, start=0, step=1):
        """
           This Hook removes the centre of mass momentum after each time step.
           **Optional arguments:**
           start
                The first iteration at which this hook should be called.
           step
                The hook will be called every `step` iterations.
        """
        self.econs_correction = 0.0
        Hook.__init__(self, start, step)

    def __call__(self, iterative):
        pass

    def init(self, iterative):
        pass

    def pre(self, iterative):
        pass

    def post(self, iterative):
        remove_com_moment(iterative.vel, iterative.masses)


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
    ff = swap_noncovalent_lammps(ff, fn_system='lammps.dat', fn_log='none', suffix='', fn_table='rcut_12.0.table', overwrite_table=False, comm=comm, move_central_cell=True)
    # Get the indexes of the ghost atoms (M atoms of TIP4P water)
    ghost_indexes = np.array([iatom for iatom in range(system.natom) if system.get_ffatype(iatom)=='TM'])

    # US MD run
    f = h5py.File('traj_1.h5', mode='w')
    hdf = HDF5Writer(f, step=100)
    #xyz = XYZWriter('traj.xyz', step=100)

    g = h5py.File('restart_1.h5', mode='w')
    hdf_restart = RestartWriter(g, step=50000)

    temp = temp*kelvin
    timestep = 0.5*femtosecond
    #press = 1e5*pascal

    thermo = NHCThermostat(temp, timecon=100.0*femtosecond)
    #baro = MTKBarostat(ff, temp, press, timecon=1000.0*femtosecond, vol_constraint=False)
    #TBC = TBCombination(thermo, baro)

    # Define CV
    N_MOF = 2208
    N_pore1 = DUMMY_NPORE1
    N_pore2 = DUMMY_NPORE2
    indices_H2O = N_MOF + 4*N_pore1 + np.array([0, 1, 2])
    indices_6MR = np.array([121, 125, 129, 135, 139, 143])
    CV = CVOopDist(ff.system, indices_H2O, indices_6MR, indices_6MR)

    # Add bias
    K = DUMMY_K*kjmol/angstrom**2
    CV0 = DUMMY_CV0*angstrom
    set_init_config_US(ff, CV0, indices_H2O)

    part_bias = ForcePartBias(ff.system)
    bias = HarmonicBias(K, CV0, CV)
    ff.add_part(part_bias)
    part_bias.add_term(bias) 

    state = [CVStateItem([CV])] #, BiasStateItem(part_bias)]

    # Add constraining sphere
    indices_pore1 = np.arange(N_MOF, N_MOF + 4*N_pore1)
    indices_pore2 = np.arange(N_MOF+4*(N_pore1+1), ff.system.natom)
    if CV0 <= -3*angstrom: indices_pore1 = np.concatenate((indices_pore1, indices_H2O))
    elif CV0 >= 3*angstrom: indices_pore2 = np.concatenate((indices_pore2, indices_H2O))

    if len(indices_pore1) > 0:
        r0 = np.dot(np.array([0.25, 0.25, 0.25]), system.cell.rvecs)
        sphere = ForcePartUpperSphere(ff.system, indices_pore1, 100*kjmol/angstrom**2, r0, 9.0*angstrom)
        ff.add_part(sphere)
    
    if len(indices_pore2) > 0:    
        r0 = np.dot(np.array([0.5, 0.5, 0.5]), system.cell.rvecs)
        sphere = ForcePartUpperSphere(ff.system, indices_pore2, 100*kjmol/angstrom**2, r0, 9.0*angstrom)
        ff.add_part(sphere)

    #ff_lammps = swap_noncovalent_lammps(ff, fn_system='lammps.dat', fn_log="log.lammps", suffix='', fn_table='rcut_12.0.table', overwrite_table=False, comm=comm)
    ghostff = GhostForceField(ff, ghost_indexes, write_ghost_pos, write_ghost_gpos) 

    vsl = VerletScreenLog(step=1000)
    md = VerletIntegrator(ghostff, timestep, state=state, hooks=[hdf, thermo, vsl, hdf_restart, RemoveComMoment()])
    md.run(1500000)

