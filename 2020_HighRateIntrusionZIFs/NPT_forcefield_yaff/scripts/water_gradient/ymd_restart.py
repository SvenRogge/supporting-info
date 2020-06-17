#!/usr/bin/env python

import os, sys
import numpy as np
import h5py as h5

from yaff import System, ForceField, HDF5Writer, RestartWriter, XYZWriter,\
    NHCThermostat, MTKBarostat, TBCombination, VerletScreenLog, VerletIntegrator, log,\
    PairPotLJ, PairPotMM3, ForcePartPair
from molmod.units import *
mpa = 1e6*pascal

from ghostatoms import GhostForceField
from tailcorr import ForcePartTailCorrection

temperatures = np.array([298.0])*kelvin
pressures = np.array([0.1,1.0])*mpa

def get_ff(rcut=12.0*angstrom, guest='water', dotail=True):
    # Load initial system
    system = System.from_file('init.chk')
    # Construct force field
    ff_fns = ['pars.txt']
    ff = ForceField.generate(system, ff_fns, rcut=rcut, alpha_scale=3.2, gcut_scale=1.5)
    # Add tail corrections
    if dotail:
        newparts = []
        for part in ff.parts:
            if part.name=='pair_mm3':
                sigmas = part.pair_pot.sigmas.copy()
                epsilons = part.pair_pot.epsilons.copy()
                onlypaulis = part.pair_pot.onlypaulis.copy()
                tr = part.pair_pot.get_truncation()
                pair_pot = PairPotMM3(sigmas, epsilons, onlypaulis, rcut, tr=tr)
                part = ForcePartPair(ff.system, ff.nlist, part.scalings, pair_pot)
                part_tail = ForcePartTailCorrection(ff.system, part.pair_pot)
                newparts.append(part_tail)
            elif part.name=='pair_lj':
                sigmas = part.pair_pot.sigmas.copy()
                epsilons = part.pair_pot.epsilons.copy()
                tr = part.pair_pot.get_truncation()
                pair_pot = PairPotLJ(sigmas, epsilons, rcut, tr=tr)
                part = ForcePartPair(ff.system, ff.nlist, part.scalings, pair_pot)
                part_tail = ForcePartTailCorrection(ff.system, part.pair_pot)
                newparts.append(part_tail)
            newparts.append(part)
        ff = ForceField(system, newparts, ff.nlist)
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
    # General settings
    rcut = 12.0*angstrom
    # Construct PES
    ghostff = get_ff(rcut=rcut)
    # Prepare MD run
    f_res = h5.File('restart_19.h5', mode='r')

    f = h5.File('traj_20.h5', mode='w')
    hdf = HDF5Writer(f, step=1000)
    g = h5.File('restart_20.h5',mode='w')
    rw = RestartWriter(g, step=100000)

    vsl = VerletScreenLog(step=5000)
    log.set_level(log.medium)
    hooks = [hdf,vsl,rw]
    md = VerletIntegrator(ghostff, hooks=hooks, restart_h5 = f_res)
    md.run(300000)
