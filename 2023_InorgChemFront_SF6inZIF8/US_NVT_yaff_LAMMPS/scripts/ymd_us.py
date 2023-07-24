#!/usr/bin/env python

import os, sys, h5py
import numpy as np
from mpi4py import MPI
from glob import glob
from yaff import *
from molmod.units import *
from molmod.constants import *


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
        diffs = np.array([self.system.pos[self.indices_plane[(i+2)%6],:]-self.system.pos[self.indices_plane[i%6],:] for i in range(len(self.indices_plane))]) # array of r13, r24, r35, r46, r51, r62, where rij = r_j - r_i. Note: takes on form 6 x n_time x 3
        normals = np.array([np.cross(diffs[i,:],diffs[(i+2)%6,:]) for i in range(len(self.indices_plane))]) # construct array of normals r13 x r35, r24 x r46, etc.
        n = np.mean(normals, axis=0)
        n = n/ np.linalg.norm(n)

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


def set_init_config_US(ff, CV0, indices, indices_plane, indices_ref):

    diffs = np.array([ff.system.pos[indices_plane[(i+2)%6],:]-ff.system.pos[indices_plane[i%6],:] for i in range(len(indices_plane))]) # array of r13, r24, r35, r46, r51, r62, where rij = r_j - r_i. Note: takes on form 6 x n_time x 3
    normals = np.array([np.cross(diffs[i,:],diffs[(i+2)%6,:]) for i in range(len(indices_plane))]) # construct array of normals r13 x r35, r24 x r46, etc.
    n = np.mean(normals, axis=0)
    n = n / np.linalg.norm(n)

    pos_SF6 = np.mean(ff.system.pos[indices], axis=0)
    pos_6MR = np.mean(ff.system.pos[indices_ref], axis=0)
    CV_old = np.dot(pos_SF6-pos_6MR,n)
    
    pos = ff.system.pos.copy()
    pos[indices] += (CV0-CV_old)*n
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

    #test_sphere()

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
    ff = swap_noncovalent_lammps(ff, fn_system='lammps.dat', fn_log='none', suffix='', fn_table='lammps_smoothei2.table', comm=comm)

    # US MD run
    f = h5py.File('traj_1.h5', mode='w')
    hdf = HDF5Writer(f, step=1000)

    g = h5py.File('restart_1.h5', mode='w')
    hdf_restart = RestartWriter(g, step=500000)

    temp = 300*kelvin
    timestep = 0.5*femtosecond
    thermo = NHCThermostat(temp, timecon=100.0*femtosecond)

    # Define CV
    indices_SF6 = np.array([2208,2209,2210,2211,2212,2213,2214])
    indices_6MR = np.array([494, 1326, 1324, 778, 768, 500])
    CV = CVOopDist(ff.system, indices_SF6, indices_6MR, indices_6MR)

    # Add bias
    K = K_DUMMY
    CV0 = CV0_DUMMY
    set_init_config_US(ff, CV0, indices_SF6, indices_6MR, indices_6MR)

    part_bias = ForcePartBias(ff.system)
    bias = HarmonicBias(K, CV0, CV)
    ff.add_part(part_bias)
    part_bias.add_term(bias) 

    state = [CVStateItem([CV])]

    vsl = VerletScreenLog(step=10000)
    md = VerletIntegrator(ff, timestep, state=state, hooks=[hdf, thermo, vsl, hdf_restart, RemoveComMoment()]) 
    md.run(2000000)
