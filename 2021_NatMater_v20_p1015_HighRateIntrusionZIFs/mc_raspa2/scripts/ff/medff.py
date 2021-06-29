#!/usr/bin/env python

import numpy as np
import h5py as h5

from yaff import ForceField, System, Scalings, Switch3, PairPotEI,\
    PairPotEiSlater1s1sCorr, ForcePartPair, NeighborList, PairPotOlpSlater1s1s, PairPotDampDisp,\
    PairPotChargeTransferSlater1s1s,\
    ForcePartEwaldReciprocal, ForcePartEwaldCorrection, ForcePartEwaldNeutralizing,\
    PairPotLJ
from yaff.pes.ext import nlist_recompute
from yaff import log
from molmod.units import angstrom, kjmol, kcalmol
from molmod.periodic import periodic

c6_table = {
    1: 6.5,
    2: 1.46,
    5: 99.5,
    6: 46.6,
    7: 24.2,
    8: 15.6,
    9: 9.5,
    10: 6.4,
    13: 528,
    14: 305.0,
    16: 134.0,
    17: 94.6,
    18: 64.2,
    30: 284.0,
    35: 162.0,
    36: 130.0,
    40: 1360.0, #10.1021/acs.jctc.6b00361
    53: np.nan,
}

alpha_table = {
    1: 4.5,
    2: 1.388,
    5: 21.0,
    6: 12.0,
    7: 7.4,
    8: 5.4,
    9: 3.8,
    10: 2.663,
    13: 60.0,
    14: 37.0,
    16: 19.6,
    17: 15.0,
    18: 11.1,
    30: 40.0,
    35: 20.0,
    36: 16.7,
    40: 112.0, #10.1021/acs.jctc.6b00361
    53: np.nan,
}



def prepare_medff(fn, atypes=None):
    '''
    Read MBIS parameters and geometry from hdf5 file
    '''
    # Load MBIS parameters
    aim_pars = {}
    print fn
    with h5.File(fn,'r') as fh5:
        for key in 'valence_charges','valence_widths','core_charges','exp_r2','exp_r4','polar_scales':
            aim_pars[key] = fh5[key][:][:]
        coords = fh5['coordinates'][:]
        numbers = fh5['numbers'][:]
        if 'rvecs' in fh5.keys():
            rvecs = fh5['rvecs'][:]
        else: rvecs = None
    # Give the system integer charge
    charges = aim_pars['valence_charges'] + aim_pars['core_charges']
    corr = np.sum(charges) - np.round( np.sum(charges) )
    aim_pars['valence_charges'][:] -= corr/numbers.shape[0]
    assert np.abs(corr) < 1e-2
    charges = aim_pars['valence_charges'] + aim_pars['core_charges']
    # Detect atom types
    if atypes is None: ffatypes = ["%s_%05d" % (periodic[numbers[iatom]].symbol,iatom) for iatom in xrange(numbers.shape[0])]
    else: ffatypes = None
    system = System(numbers, coords, rvecs=rvecs, ffatypes=ffatypes)
    system.detect_bonds()
    if atypes is not None: 
        system.detect_ffatypes(atypes)
    return system, aim_pars

def get_dispersion_coefficients(valence_widths, polar_scales, exp_r2, exp_r4, numbers):
    natom = valence_widths.shape[0]
    c6s = np.zeros((natom, natom))
    c8s = np.zeros((natom, natom))
    alphas = np.zeros((natom,))
    for iatom0 in xrange(natom):
        alpha0 = alpha_table[numbers[iatom0]]*polar_scales[iatom0]
        alphas[iatom0] = alpha0
        c60 = c6_table[numbers[iatom0]]*polar_scales[iatom0]**2
        for iatom1 in xrange(iatom0+1):
            alpha1 = alpha_table[numbers[iatom1]]*polar_scales[iatom1]
            c61 = c6_table[numbers[iatom1]]*polar_scales[iatom1]**2
            c6 = 2.0*c60*c61/(alpha1/alpha0*c60 + alpha0/alpha1*c61)
            c8 = 1.5*(exp_r4[iatom0]/exp_r2[iatom0] + exp_r4[iatom1]/exp_r2[iatom1])*c6
            c6s[iatom0,iatom1] = c6
            c6s[iatom1,iatom0] = c6
            c8s[iatom0,iatom1] = c8
            c8s[iatom1,iatom0] = c8
    R_cross = np.tile(valence_widths,natom).reshape((natom,natom))
    R_cross = 2.0/(R_cross + R_cross.T)
    return c6s, c8s, R_cross



neigh_dtype = [
    ('a', int), ('b', int), ('d', float),        # a & b are atom indexes, d is the distance
    ('dx', float), ('dy', float), ('dz', float), # relative vector (includes cell vectors of image cell)
    ('r0', int), ('r1', int), ('r2', int)        # position of image cell.
]

class ForcePartEwaldReciprocalGuest(ForcePartEwaldReciprocal):
    def prepare(self, nguest=0):
        self.nguest = nguest
        self.cosfacs = np.zeros(2*self.gmax+1)
        self.sinfacs = np.zeros(2*self.gmax+1)
        self.prefacs = np.zeros(2*self.gmax+1)
        self.kvectors = np.zeros((2*self.gmax[0]+1, 2*self.gmax[1]+1, 2*self.gmax[2]+1, 3))
        self.dots = np.zeros((2*self.gmax[0]+1, 2*self.gmax[1]+1, 2*self.gmax[2]+1, self.nguest+1))
        self.cosfac_guests = np.zeros(2*self.gmax+1)
        self.sinfac_guests = np.zeros(2*self.gmax+1)
        self.kvecs = self.system.cell.gvecs*2.0*np.pi
        fac1 = 8.0*np.pi/self.system.cell.volume
        fac2 = 0.25/self.alpha/self.alpha
        kcut = (self.gcut*2.0*np.pi)**2
        energy = 0.0
        for g0 in xrange(-self.gmax[0],self.gmax[0]+1):
            for g1 in xrange(-self.gmax[1],self.gmax[1]+1):
                for g2 in xrange(0,self.gmax[2]+1):
                    if g2==0:
                        if g1<0: continue
                        if g1==0 and g0<=0: continue
                    k = np.dot(self.kvecs, np.array([g0,g1,g2]))
                    ksq = np.dot(k,k)
                    if ksq > kcut: continue
                    index = self.gmax[0]+g0,self.gmax[1]+g1,self.gmax[2]+g2
                    self.kvectors[index] = k
                    self.prefacs[index] = np.exp(-fac2*ksq)/ksq*fac1
#                    for iatom in xrange(self.nguest+1,self.system.natom):
#                        self.cosfacs[index] += self.system.charges[iatom]*np.cos(np.dot(k,self.system.pos[iatom])) 
#                        self.sinfacs[index] += self.system.charges[iatom]*np.sin(np.dot(k,self.system.pos[iatom])) 
#                    self.cosfacs[index] *= np.exp(-fac2*ksq)/ksq*fac1
#                    self.sinfacs[index] *= np.exp(-fac2*ksq)/ksq*fac1
        dots = np.einsum('ghja,ia->ghji',self.kvectors, self.system.pos[self.nguest+1:])
        self.cosfacs = np.einsum('i,ghji',self.system.charges[self.nguest+1:],np.cos(dots))
        self.sinfacs = np.einsum('i,ghji',self.system.charges[self.nguest+1:],np.sin(dots))
        self.cosfacs *= self.prefacs
        self.sinfacs *= self.prefacs


    def _internal_compute(self, gpos, vtens, hess):
        self.dots[:] = np.einsum('ghja,ia->ghji',self.kvectors, self.system.pos[:self.nguest+1])
        self.cosfac_guests[:] = np.einsum('i,ghji',self.system.charges[:self.nguest+1],np.cos(self.dots))
        self.sinfac_guests[:] = np.einsum('i,ghji',self.system.charges[:self.nguest+1],np.sin(self.dots))
        return np.einsum('ghj,ghj',self.cosfac_guests,self.cosfacs) + np.einsum('ghj,ghj',self.sinfac_guests,self.sinfacs)

class XmerNeighborList(NeighborList):
    '''
    Neighborlist where only interactions between molecules are considered.
    This means the list never has to be rebuilt, only recomputed
    '''
    def __init__(self, system, group_ids):
        # Initialize
        self.system = system
        self.skin = 0
        self.rcut = 0.0
        self.rmax = None
        self._pos_old = system.pos.copy()
        # Construct neighbor list once and for all
        self.nneigh = system.natom**2
        for group in set(group_ids):
            self.nneigh -= np.sum(group_ids==group)**2
        self.nneigh /= 2
        self.neighs = np.empty(self.nneigh, dtype=neigh_dtype)
        index = 0
        for i in xrange(system.natom):
            for j in xrange(i):
                if group_ids[i]!=group_ids[j]:
                    self.neighs[index]['a'] = i
                    self.neighs[index]['b'] = j
                    for k in xrange(3): self.neighs[index]['r%d'%k] = 0
                    index += 1
        self.update()
        self.nguest = np.sum(group_ids==0)-1

    def update(self):
        nlist_recompute(self.system.pos, self._pos_old, self.system.cell, self.neighs)

    def _need_rebuild(self):
        return False



class MEDFFBase(object):
    def __init__(self, system, aim_data, upars, covalent=None, olp_induction=False):
        self.system = system
        self.upars = upars
        self.covalent = covalent
        self.olp_induction = olp_induction
        # Add aim data to system
        self.system.charges = aim_data['core_charges'] + aim_data['valence_charges']
        self.system.slater1s_N = aim_data['valence_charges']
        self.system.slater1s_Z = aim_data['core_charges']
        self.system.slater1s_widths = aim_data['valence_widths']
        self.components = ['medff-electrostatics','medff-overlap','medff-disp6','medff-disp8']
#        self.nai = -aim_data['valence_charges']-aim_data['valences']

    def update_pos(self, pos):
        self.ff.update_pos(pos)

    def compute(self, gpos=None, vtens=None):
        #self.system.pos[:] = pos
        self.nlist.update()
        self.ff.compute(gpos, vtens)
        return self.ff.energy
        energies = np.zeros((len(self.components)))
        for part in self.ff.parts:
            if part.name.startswith('pair_ei') or part.name.startswith('ewald'):
                index = self.components.index('medff-electrostatics')
            elif part.name == "pair_olpslater1s1s":
                index = self.components.index('medff-overlap')
            elif part.name == "pair_dampdisp_6":
                index = self.components.index('medff-disp6')
            elif part.name == "pair_dampdisp_8":
                index = self.components.index('medff-disp8')
            else: raise NotImplementedError
            energies[index] += part.energy
#            print "Adding energy of %25s (%+8.1e) to %15s" % (part.name, part.energy, self.components[index])
        return energies
#        assert False
        return [part.energy for part in self.ff.parts]

    def construct_ff(self, upars, aim_data):
        # Add all terms to the force field
        self.parts = []
        self._init_electrostatics()     
        self._init_dispersion()
        self._init_exchange()
        if not self.olp_induction:
            self._init_induction()
        if self.covalent is not None:
            # Load bonded forcefield from QuickFF
            ff_valence = ForceField.generate(self.system,self.covalent)
            self.parts.append(ff_valence.part_valence)
        # Construct Yaff force field
        self.ff = ForceField(self.system, self.parts, self.nlist)

    def _init_electrostatics(self):
        # Point charges
        pair_pot_ei_real  = PairPotEI(self.system.charges, 0.0, self.rcut, tr=self.truncation)
        part_pair_ei_real = ForcePartPair(self.system, self.nlist, self.scalings, pair_pot_ei_real)
        self.parts.append(part_pair_ei_real)
        # Add pair potential for Slater corrections
        pair_pot_ei_slater = PairPotEiSlater1s1sCorr(self.system.slater1s_widths, self.system.slater1s_N, self.system.slater1s_Z, self.rcut, tr=self.truncation)
        part_pair_ei_slater = ForcePartPair(self.system, self.nlist, self.scalings, pair_pot_ei_slater)
        self.parts.append(part_pair_ei_slater)

    def _init_exchange(self):
        upar = self.upars['ex_scale']
        if self.olp_induction: upar -= self.upars['ct_scale']
        pair_pot_ex  = PairPotOlpSlater1s1s(self.system.slater1s_widths, self.system.slater1s_N, upar, self.rcut, tr=self.truncation)
        part_pair_ex = ForcePartPair(self.system, self.nlist, self.scalings, pair_pot_ex)
        self.parts.append(part_pair_ex)

    def _init_induction(self):
        self.apars = np.zeros((self.system.natom), float)
        self.bpars = np.zeros((self.system.natom), float)
        for iatom in xrange(self.system.natom):
            Zi = self.system.numbers[iatom]
            self.apars[iatom] = self.upars['A_%s'%(periodic[Zi].symbol)]
            self.bpars[iatom] = self.upars['B_%s'%(periodic[Zi].symbol)]
        pair_pot_ct = PairPotChargeTransferSlater1s1s(self.system.slater1s_widths, self.system.slater1s_N, self.system.charges, self.apars, self.bpars, self.upars['Uind'], self.rcut, tr=self.truncation)
        part_pair_ct = ForcePartPair(self.system, self.nlist, self.scalings, pair_pot_ct)
        self.parts.append(part_pair_ct) 

    def get_dispersion_coefficients(self, scale_c6, scale_c8, valence_widths, polar_scales, exp_r2, exp_r4, numbers):
        natom = valence_widths.shape[0]
        #ffatypes_disp = np.arange(natom)
        c6s = np.zeros((natom, natom))
        c8s = np.zeros((natom, natom))
        self.alphas = np.zeros((natom,))
        for iatom0 in xrange(natom):
            alpha0 = alpha_table[numbers[iatom0]]*polar_scales[iatom0]
            self.alphas[iatom0] = alpha0
            c60 = c6_table[numbers[iatom0]]*polar_scales[iatom0]**2
            for iatom1 in xrange(iatom0+1):
                alpha1 = alpha_table[numbers[iatom1]]*polar_scales[iatom1]
                c61 = c6_table[numbers[iatom1]]*polar_scales[iatom1]**2
                c6 = 2.0*c60*c61/(alpha1/alpha0*c60 + alpha0/alpha1*c61)
                c8 = 1.5*(exp_r4[iatom0]/exp_r2[iatom0] + exp_r4[iatom1]/exp_r2[iatom1])*c6
                c6s[iatom0,iatom1] = c6
                c6s[iatom1,iatom0] = c6
                c8s[iatom0,iatom1] = c8
                c8s[iatom1,iatom0] = c8
        c6s *= scale_c6
        c8s *= scale_c8
        R_cross = np.tile(valence_widths,natom).reshape((natom,natom))
        R_cross = 2.0/(R_cross + R_cross.T)
        return c6s, c8s, R_cross

    def _init_dispersion(self):
        pair_pot_disp6 = PairPotDampDisp(self.ffatypes_disp,self.c6s,self.R_cross,self.rcut,tr=self.truncation, power=6)
        part_pair_disp6 = ForcePartPair(self.system, self.nlist, self.scalings, pair_pot_disp6)
        part_pair_disp6.name += '_6'
        self.parts.append(part_pair_disp6)
        pair_pot_disp8 = PairPotDampDisp(self.ffatypes_disp,self.c8s,self.R_cross,self.rcut,tr=self.truncation, power=8)
        part_pair_disp8 = ForcePartPair(self.system, self.nlist, self.scalings, pair_pot_disp8)
        part_pair_disp8.name += '_8'
        self.parts.append(part_pair_disp8)

class MEDFFLiquid(MEDFFBase):
    def __init__(self, system, aim_data, upars, covalent=None, olp_induction=False, scales=[0.0,0.0,0.0],rcut=15.0*angstrom,nguest=-1, do_trunc=True):
        self.nmol = system.natom/aim_data['core_charges'].shape[0]
        self.nguest = nguest
#        print "NMOL", self.nmol
        assert system.natom == aim_data['core_charges'].shape[0]*self.nmol
        self.monoatom = system.natom/self.nmol
        for prop in 'core_charges','valence_charges','valence_widths':
            aim_data[prop] = np.tile(aim_data[prop], self.nmol)
        MEDFFBase.__init__(self, system, aim_data, upars, covalent=covalent, olp_induction=olp_induction)
        self.rcut = rcut
        self.scalings = Scalings(system, scales[0], scales[1], scales[2])
        if do_trunc:
            self.truncation = Switch3(1.4*angstrom)
        else:
#            print "Not applying truncation"
            self.truncation = None
        self.nlist = NeighborList(self.system, 0.0, nguest=nguest)
        self.ffatypes_disp = np.tile(np.arange(self.monoatom),self.nmol)
        self.c6s, self.c8s, self.R_cross = self.get_dispersion_coefficients(upars['scale_c6'], upars['scale_c8'],aim_data['valence_widths'][:self.monoatom],
                aim_data['polar_scales'],aim_data['exp_r2'],aim_data['exp_r4'], self.system.numbers[:self.monoatom])
        self.construct_ff(upars, aim_data)

    def _init_electrostatics(self):
        # Point charges
        assert self.system.cell.nvec == 3
        alpha = 3.5/self.rcut
        if self.nguest>-1:
            part_ewald_reci = ForcePartEwaldReciprocalGuest(self.system, alpha, gcut=2.5*alpha)
            part_ewald_reci.prepare(self.nguest)
            self.parts.append(part_ewald_reci)
        else:
            part_ewald_reci = ForcePartEwaldReciprocal(self.system, alpha, gcut=2.5*alpha)
            part_ewald_corr = ForcePartEwaldCorrection(self.system, alpha, self.scalings)
            part_ewald_neut = ForcePartEwaldNeutralizing(self.system, alpha)
            self.parts.append(part_ewald_reci)
            self.parts.append(part_ewald_corr)
            self.parts.append(part_ewald_neut)
        pair_pot_ei_real  = PairPotEI(self.system.charges, alpha, self.rcut, tr=self.truncation)
        part_pair_ei_real = ForcePartPair(self.system, self.nlist, self.scalings, pair_pot_ei_real)
        self.parts.append(part_pair_ei_real)
        # Add pair potential for Slater corrections
        pair_pot_ei_slater = PairPotEiSlater1s1sCorr(self.system.slater1s_widths, self.system.slater1s_N, self.system.slater1s_Z, self.rcut, tr=self.truncation)
        part_pair_ei_slater = ForcePartPair(self.system, self.nlist, self.scalings, pair_pot_ei_slater)
        self.parts.append(part_pair_ei_slater)

class MEDFFMonomer(MEDFFBase):
    def __init__(self, system, aim_data, upars, covalent=None, olp_induction=False, scales=[0.0,0.0,0.0], rcut=200.0):
        self.group_ids = np.zeros(system.natom,dtype=int)
        MEDFFBase.__init__(self, system, aim_data, upars, covalent=covalent, olp_induction=olp_induction)
        # Force-field parameters
        self.rcut = rcut
        self.scalings = Scalings(system, scales[0], scales[1], scales[2])
        self.truncation = Switch3(2.0*angstrom)
        self.nlist = NeighborList(self.system)
        self.ffatypes_disp = np.arange(self.system.natom)
        self.c6s, self.c8s, self.R_cross = self.get_dispersion_coefficients(upars['scale_c6'], upars['scale_c8'],aim_data['valence_widths'],
                aim_data['polar_scales'],aim_data['exp_r2'],aim_data['exp_r4'],self.system.numbers)
        self.construct_ff(upars, aim_data)

class MEDFFDimer(MEDFFBase):
    def __init__(self, system, group_ids, aim_data, upars, olp_induction=False, scales=[0.0,0.0,0.0], rcut=150.0*angstrom, do_trunc=True):
        self.group_ids = group_ids
        MEDFFBase.__init__(self, system, aim_data, upars, olp_induction=olp_induction)
        # Force-field parameters
        self.rcut = rcut
        self.scalings = Scalings(system, scales[0], scales[1], scales[2])
        if do_trunc:
            self.truncation = Switch3(1.4*angstrom)
        else:
            self.truncation = None
        #self.nlist = NeighborList(self.system)#, self.group_ids)
        self.nlist = XmerNeighborList(self.system, self.group_ids)
        self.ffatypes_disp = np.arange(self.system.natom)
        self.c6s, self.c8s, self.R_cross = self.get_dispersion_coefficients(upars['scale_c6'], upars['scale_c8'],aim_data['valence_widths'],
                aim_data['polar_scales'],aim_data['exp_r2'],aim_data['exp_r4'],self.system.numbers)
        self.construct_ff(upars, aim_data)

    def compute(self, gpos=None, vtens=None):
        #self.system.pos[:] = pos
        self.nlist.update()
        self.ff.compute(gpos, vtens)
        return self.ff.energy
        energies = np.zeros((len(self.components)))
        for part in self.ff.parts:
            if part.name.startswith('pair_ei') or part.name.startswith('ewald'):
                index = self.components.index('medff-electrostatics')
            elif part.name == "pair_olpslater1s1s":
                index = self.components.index('medff-overlap')
            elif part.name == "pair_dampdisp_6":
                index = self.components.index('medff-disp6')
            elif part.name == "pair_dampdisp_8":
                index = self.components.index('medff-disp8')
            else: raise NotImplementedError
            energies[index] += part.energy
        return energies
        return [part.energy for part in self.ff.parts]

    def get_long_range_coefficients(self):
        natom = self.system.natom
        mask = np.tile(self.group_ids, natom).reshape((natom,natom))
        mask = mask+mask.T==1
        # Dispersion
        C6 = -np.sum(self.ff.parts[2].pair_pot.cn_cross[mask])/2.0
        C8 = -np.sum(self.ff.parts[3].pair_pot.cn_cross[mask])/2.0
        return [[6,C6], [8,C8]]
