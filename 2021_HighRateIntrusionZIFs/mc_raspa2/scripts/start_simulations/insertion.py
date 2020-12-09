#!/usr/bin/env python

import os
import numpy as np
import h5py as h5
#np.random.seed(34)

from yaff import System, ForceField, log
log.set_level(log.silent)
from molmod.units import angstrom, kjmol, bar, kelvin, liter, centimeter, gram, meter, kilogram, atm
from molmod.constants import boltzmann, avogadro

from test_ghostatoms import GhostForceField, write_ghost_positions_n2

ffname = 'qff-dreiding-mbis-tip4p'
fw = 'zif8-jelle'
rcut = 12.0*angstrom
ff_fns = ['../ffpars/%s/%s/pars.txt'%(ffname,fw)]

def rand_rotation_matrix(deflection=1.0, randnums=None):
    """
    Creates a random rotation matrix.
    
    deflection: the magnitude of the rotation. For 0, no rotation; for 1, competely random
    rotation. Small deflection => small perturbation.
    randnums: 3 random numbers in the range [0, 1]. If `None`, they will be auto-generated.
    """
    # from http://www.realtimerendering.com/resources/GraphicsGems/gemsiii/rand_rotation.c
    
    if randnums is None:
        randnums = np.random.uniform(size=(3,))
        
    theta, phi, z = randnums
    
    theta = theta * 2.0*deflection*np.pi  # Rotation about the pole (Z).
    phi = phi * 2.0*np.pi  # For direction of pole deflection.
    z = z * 2.0*deflection  # For magnitude of pole deflection.
    
    # Compute a vector V used for distributing points over the sphere
    # via the reflection I - V Transpose(V).  This formulation of V
    # will guarantee that if x[1] and x[2] are uniformly distributed,
    # the reflected points will be uniform on the sphere.  Note that V
    # has length sqrt(2) to eliminate the 2 in the Householder matrix.
    
    r = np.sqrt(z)
    Vx, Vy, Vz = V = (
        np.sin(phi) * r,
        np.cos(phi) * r,
        np.sqrt(2.0 - z)
        )
    
    st = np.sin(theta)
    ct = np.cos(theta)
    
    R = np.array(((ct, st, 0), (-st, ct, 0), (0, 0, 1)))
    
    # Construct the rotation matrix  ( V Transpose(V) - I ) R.
    
    M = (np.outer(V, V) - np.eye(3)).dot(R)
    return M

def gen_water():
    pos = np.array([[0.00000000,0.00000000,0.00000000],
                    [0.7570,0.5859,0.00000000],
                    [-0.7570,0.5859,0.00000000],
                    [0.00000000,0.1546,0.00000000]])*angstrom
    M = rand_rotation_matrix()
    return np.einsum('ab,ib',M,pos).T

def add_molecule(system, guest):
    if guest=='n2':
        natom = 3
    elif guest=='ar':
        natom = 1
    elif guest=='water':
        natom = 4
    else: raise NotImplementedError
    fwcharges = system.charges.copy()
    ffold = ForceField.generate(system, ff_fns, rcut=rcut)
    ffold.system.charges[:] = fwcharges
    eold = ffold.compute()
    newnatom = system.natom + natom
    newpos = np.zeros((newnatom,3))
    newpos[:system.natom] = system.pos
    newnumbers = np.zeros((newnatom), dtype=int)
    newffatypes = list(system.ffatypes[system.ffatype_ids])
    newnumbers[:system.natom] = system.numbers
    newbonds = list(system.bonds)
    if guest=='n2':
        newnumbers[system.natom:] = [7,7,99]
        for ffa in ['NN2','NN2','GN2']:
            newffatypes.append(ffa)
        for new in ([0,1],[0,2],[1,2]):
            newbonds.append([new[0]+system.natom,new[1]+system.natom])
    elif guest=='ar':
        newnumbers[system.natom:] = [18]
        newffatypes.append('Ar')
    elif guest=='water':
        newnumbers[system.natom:] = [8,1,1,98]
        for ffa in ['TO','TH','TH','TM']:
            newffatypes.append(ffa)
        monobonds = [[0,1],[0,2],[0,3]]
        for new in monobonds:
            newbonds.append([new[0]+system.natom,new[1]+system.natom])
        newcharges = np.array([0.0,0.5564,0.5564,-2.0*0.5564])
    else: raise NotImplementedError
    newbonds = np.array(newbonds, dtype=int)
    newsystem = System(newnumbers, newpos, ffatypes=newffatypes, rvecs=system.cell.rvecs, bonds=newbonds)
    ff = ForceField.generate(newsystem, ff_fns, rcut=rcut)
    ff.system.charges[:fwcharges.shape[0]] = fwcharges
    ff.system.charges[fwcharges.shape[0]:] = newcharges

    monosystem = System(newnumbers[system.natom:], newpos[system.natom:], ffatypes=newffatypes[system.natom:], rvecs=system.cell.rvecs, bonds=monobonds)
    monoff = ForceField.generate(monosystem, ff_fns, rcut=rcut)
    monoff.system.charges[:] = newcharges

    imax = 10000
    for counter in xrange(imax):
        pos = np.random.rand(3)
        pos = np.dot(system.cell.rvecs.T, pos)
        newsystem.pos[system.natom] = pos
        if guest=='n2':
            # WARNING! Random rotations not uniformly sampled here!
            # only use for stuff like initial configurations
            u = np.random.rand(3)
            u /= np.linalg.norm(u)
            newsystem.pos[system.natom+1] = pos + u*1.091163*angstrom
            newsystem.pos[system.natom+2] = pos + 0.5*u*1.091163*angstrom
        elif guest=='ar': pass
        elif guest=='water':
            newsystem.pos[system.natom:] = pos + gen_water()
        else: raise NotImplementedError
        ff.update_pos(newsystem.pos)
        enew = ff.compute()
        monoff.update_pos(newsystem.pos[system.natom:])
        enew -= monoff.compute()
        if enew-eold < 0.0*kjmol:
            print "Attempt %5d Delta = %12.3f kJ/mol" % (counter,(enew-eold)/kjmol)
            break
        if counter==imax-1:
            assert False
    return newsystem

if __name__=='__main__':
    guest = 'water'
    for nmol in xrange(4,81,4):
        system = System.from_file('../ffpars/%s/%s/optcell.chk'%(ffname,fw))
        workdir = '../insertion/%s/%s/%s/%05d'% (fw,ffname,guest, nmol)
        if not os.path.isdir(workdir): os.makedirs(workdir)
        mass_fw = np.sum(system.masses)
        rho = mass_fw/system.cell.volume
        loading = nmol*boltzmann*273.15*kelvin/atm/mass_fw
        print "="*100
        print "Inserting %5d %s molecules corresponding to loading of %8.2f cm**3(STP)/g = %8.2f mmol/g" % (nmol,guest,loading/centimeter**3*gram,nmol/mass_fw/avogadro*gram*1000.0)
        for imol in xrange(nmol):
            system = add_molecule(system, guest)
        system.to_file(os.path.join(workdir,'init.chk'))
