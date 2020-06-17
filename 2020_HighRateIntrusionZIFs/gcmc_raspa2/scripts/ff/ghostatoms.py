#!/usr/bin/env python

import numpy as np
np.random.seed(305)

from yaff import System, ForceField, NeighborList, Scalings, PairPotLJ,\
 ForcePartPair, ForcePart, log
from molmod.units import angstrom, kcalmol

class GhostForceField(ForceField):
    '''A force field that implicitly contains ghost atoms. Ghost atoms do not
    play a role in sampling (verlet, optimization, ) but influence the energy.
    The position of ghost atoms is based on geometric rules.'''
    def __init__(self, ff, ghost_indexes, write_ghost_pos, write_ghost_gpos):
        """

           **Arguments**

           ff
                An instance of the ``ForceField`` class, containing all atoms
                (including ghost atoms). This is used for the energy evaluation

           ghost_indexes
                NumPy array containing the indexes of the ghost atoms

           write_ghost_pos
                A function that updates the positions of the ghost atoms based
                on the positions of the other atoms. The arguments of this
                function are a NumPy array with positions of ALL atoms and the
                indexes of the ghost atoms. This function writes directly to
                the NumPy array with positions.

        """
        self.ghost_indexes = ghost_indexes
        self.atom_indexes = np.array([iatom for iatom in xrange(ff.system.natom) if not iatom in self.ghost_indexes], dtype=int)
        self.system = ff.system.subsystem(self.atom_indexes)
        self.write_ghost_pos = write_ghost_pos
        self.write_ghost_gpos = write_ghost_gpos
        ForcePart.__init__(self, 'ghostff', self.system)
        self.ff = ff
        self.gpos_full = np.zeros((self.ff.system.pos.shape[0], self.ff.system.pos.shape[1]))
        self.vtens_full = np.zeros((3,3))
        self.parts = [self.ff]
        # These attributes need to be defined but are never used...
        self.gpos_parts = None
        self.vtens_parts = None
        if log.do_medium:
            with log.section('FFINIT'):
                log('Force field with %d ghost atoms and %d real atoms.' % (self.ghost_indexes.shape[0], self.atom_indexes.shape[0]))

    def update_rvecs(self, rvecs):
        '''See :meth:`yaff.pes.ff.ForcePart.update_rvecs`'''
        ForcePart.update_rvecs(self, rvecs)
        ForcePart.update_rvecs(self.ff, rvecs)
        self.system.cell.update_rvecs(rvecs)
        self.ff.system.cell.update_rvecs(rvecs)
        if self.ff.nlist is not None:
            self.ff.nlist.update_rmax()
            self.ff.needs_nlist_update = True

    def update_pos(self, pos):
        '''See :meth:`yaff.pes.ff.ForcePart.update_pos`'''
        ForcePart.update_pos(self, pos)
        ForcePart.update_pos(self.ff, pos)
        self.system.pos[:] = pos
        self.ff.system.pos[self.atom_indexes] = pos
        # Call function that computes the positions of the ghost atoms based
        # on position of other atoms
        self.write_ghost_pos(self.ff.system, self.ghost_indexes)
        if self.ff.nlist is not None:
            self.ff.needs_nlist_update = True

    def _internal_compute(self, gpos, vtens):
        if gpos is None:
            my_gpos = None
        else:
            my_gpos = self.gpos_full
            my_gpos[:] = 0.0
        if vtens is None:
            my_vtens = None
        else:
            if gpos is None: raise NotImplementedError("Cannot compute vtens without gpos")
            my_vtens = self.vtens_full
            my_vtens[:] = 0.0
        result = self.ff.compute(my_gpos, vtens)
        if not ((gpos is None) and (vtens is None)):
            if gpos is not None and np.isnan(my_gpos).any():
                raise ValueError('Some gpos element(s) is/are not-a-number (nan).')
            self.write_ghost_gpos(my_gpos, my_vtens, self.ff.system, self.ghost_indexes)
            if gpos is not None:
                gpos[:] += my_gpos[self.atom_indexes]
            if vtens is not None:
                vtens[:] += my_vtens
        return result

def get_n2_system(nmol, L):
    d = 1.45*angstrom # N2 bond length
    natom = 3   # Number of sites in each molecule
    # Give molecules random position and orientation
    pos = np.zeros((nmol*natom,3))
    bonds = []
    for imol in xrange(nmol):
        pos[imol*natom] = np.random.rand(3)*L
        # WARNING! Random rotations not uniformly sampled here!
        u = np.random.rand(3)
        u /= np.linalg.norm(u)
        pos[imol*natom+1] = pos[imol*natom] + u*d
        # Ghost atom is at COM
        pos[imol*natom+2] = 0.5*(pos[imol*natom]+pos[imol*natom+1])
        bonds.append([imol*natom,imol*natom+1])
        bonds.append([imol*natom,imol*natom+2])
        bonds.append([imol*natom+1,imol*natom+2])
    bonds = np.asarray(bonds)
    # The ghost atom gets atomic number 99, please do not use Einsteinium(99)
    # in these simulations :)
    numbers = np.tile([7,7,99], nmol)
    ffatypes = np.tile(['NN2','NN2','GN2'],nmol)
    ghost_indexes = np.arange(2,natom*nmol,natom)
    rvecs = np.eye(3)*L
    # Construct the Yaff system including ghost atoms
    system = System(numbers, pos, rvecs=rvecs, bonds=bonds, ffatypes=ffatypes)
    return system, ghost_indexes

def get_n2_forcefield(system):
    # Construct a force field for the system including ghost atoms
    nlist = NeighborList(system)
    scalings = Scalings(system, 0.0, 0.0, 1.0)
    # Initialize (random) parameters
    rminhalf_table = {
        7: 1.7000*angstrom,
        99: 1.7682*angstrom
    }
    epsilon_table = {
        7: -0.1970*kcalmol,
        99: -0.1521*kcalmol,
    }
    sigmas = np.zeros(system.natom, float)
    epsilons = np.zeros(system.natom, float)
    for i in xrange(system.natom):
        sigmas[i] = rminhalf_table[system.numbers[i]]*(2.0)**(5.0/6.0)
        epsilons[i] = epsilon_table[system.numbers[i]]
    # Construct the pair potential and part
    pair_pot = PairPotLJ(sigmas, epsilons, 15*angstrom)
    part_pair = ForcePartPair(system, nlist, scalings, pair_pot)
    ff = ForceField(system, [part_pair], nlist=nlist)
    return ff

def write_ghost_positions_n2(pos, ghost_indexes):
    # For N2, each ghost atom is at the center between the two atoms that
    # precede the ghost atom in the system
    for iatom in ghost_indexes:
        pos[iatom] = 0.5*(pos[iatom-2]+pos[iatom-1])

def write_ghost_forces_n2(gpos, ghost_indexes):
    # For N2, each ghost atom is at the center between the two atoms that
    # precede the ghost atom in the system
    for iatom in ghost_indexes:
        gpos[iatom-2] += 0.5*gpos[iatom]
        gpos[iatom-1] += 0.5*gpos[iatom]

if __name__=='__main__':
    L = 20.0*angstrom # Box size
    nmol = 2    # Number of N2 molecules
    # Construction of system and force field including ghost atoms
    system, ghost_indexes = get_n2_system(nmol, L)
    ff_full = get_n2_forcefield(system)
    # Force field where ghosts are not implicitly included, but still taken
    # into account for force calculations
    ff = GhostForceField(ff_full, ghost_indexes, write_ghost_positions_n2, write_ghost_forces_n2)
    # Tests
    gpos = np.zeros((ff.system.pos.shape[0],ff.system.pos.shape[1]))
    gpos_full = np.zeros((ff_full.system.pos.shape[0],ff_full.system.pos.shape[1]))
    e = ff.compute(gpos=gpos)
    e_full = ff_full.compute(gpos=gpos_full)
    assert e==e_full
    assert np.all(np.abs(np.sum(gpos, axis=0)) < 1e-16)
    newpos = ff.system.pos.copy()*(1.0+np.random.normal(0.0,0.05*angstrom,(ff.system.natom,3)))
    ff.update_pos(newpos)
    gpos[:] = 0.0
    gpos_full[:] = 0.0
    e = ff.compute(gpos=gpos)
    e_full = ff_full.compute(gpos=gpos_full)
    assert e==e_full
    assert np.all(np.abs(np.sum(gpos, axis=0)) < 1e-16)
    print e
