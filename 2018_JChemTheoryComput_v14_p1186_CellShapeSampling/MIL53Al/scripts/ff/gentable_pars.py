#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import glob #To use *.chk as filename

#Append lammps folder to path
sys.path.append("ff-lammps") 

from yaff import System, ForceField, log
from molmod.units import angstrom, kjmol

from mylammps import *

def write_lammps_table_jelle(ff):
    from yaff import NeighborList, PairPotEI, Switch3, ForcePartPair, Scalings,\
        PairPotMM3
    nffa = ff.system.ffatypes.shape[0]
    numbers = []
    for i in xrange(nffa):
        index0 = np.where(ff.system.ffatype_ids==i)[0][0]
        if np.sum(ff.system.ffatype_ids==i)>1:
            index1 = np.where(ff.system.ffatype_ids==i)[0][1]
            numbers.append( [index0, index1] )
        for j in xrange(i+1,nffa):
            index1 = np.where(ff.system.ffatype_ids==j)[0][0]
            numbers.append( [index0, index1] )
    part_names = [part.name for part in ff.parts]
    print part_names
    ftab = open('data-lammps/lammps_smoothei2.table','w')
    ftab.write("# LAMMPS tabulated potential generated by Yaff\n")
    ftab.write("# All quantities in atomic units\n")
    ftab.write("# The names of the tables refer to the ffatypes that have to be used in the Yaff system\n")
    ftab.write("#%4s %13s %21s %21s\n" % ("i","d","V","F"))

    for number in numbers:
        pos = np.zeros((2,3))
        pos[1,2] = 1.0
        smallsys = System(system.numbers[number],pos,charges=system.charges[number],radii=system.radii[number])
        smallsys_noradii = System(system.numbers[number],pos,charges=system.charges[number])
        tr = ff.parts[2].pair_pot.get_truncation()
        rcut = ff.parts[2].pair_pot.rcut
        alpha = ff.parts[2].pair_pot.alpha
        print "alpha:",alpha
	nlist = NeighborList(smallsys)
        pair_pot = PairPotEI(smallsys.charges, alpha, rcut, tr=Switch3(tr.width), radii=smallsys.radii)
#        pair_pot = PairPotEI(smallsys.charges, 0.0, rcut, tr=None, radii=smallsys.radii)
        scalings = Scalings(smallsys, scale1=1.0, scale2=1.0, scale3=1.0)
        part0 = ForcePartPair(smallsys,nlist,scalings,pair_pot)
        ff0 = ForceField(smallsys,[part0],nlist)
        ff0.compute()

        nlist1 = NeighborList(smallsys_noradii)
#        pair_pot = PairPotEI(smallsys_noradii.charges, alpha, rcut, tr=Switch3(tr.width))
        pair_pot = PairPotEI(smallsys_noradii.charges, alpha, rcut, tr=None)
        scalings = Scalings(smallsys_noradii, scale1=1.0, scale2=1.0, scale3=1.0)
        part1 = ForcePartPair(smallsys_noradii,nlist1,scalings,pair_pot)
        ff1 = ForceField(smallsys_noradii,[part1],nlist1)
        ff1.compute()

        tr = ff.parts[1].pair_pot.get_truncation()
        rcut = ff.parts[1].pair_pot.rcut
        nlist2 = NeighborList(smallsys)
        pair_pot = PairPotMM3(ff.parts[1].pair_pot.sigmas[number],ff.parts[1].pair_pot.epsilons[number],ff.parts[1].pair_pot.onlypaulis[number],
                rcut, Switch3(tr.width) )
        scalings = Scalings(smallsys, scale1=1.0, scale2=1.0, scale3=1.0)
        part2 = ForcePartPair(smallsys,nlist,scalings,pair_pot)
        ff2 = ForceField(smallsys,[part2],nlist)
        ff2.nlist.nneigh = 1
        ff2.compute()

        distances = np.linspace(0.5*angstrom, rcut, 5000)
        energies = []
        toplot = []
        for d in distances:
            gposnn0 = np.zeros(ff0.system.pos.shape, float)
            ff0.nlist.neighs[0] = (0, 1, d, 0.0, 0.0, d, 0, 0, 0)
            energy0 = ff0.compute(gpos=gposnn0)
            gposnn1 = np.zeros(ff1.system.pos.shape, float)
            ff1.nlist.neighs[0] = (0, 1, d, 0.0, 0.0, d, 0, 0, 0)
            energy1 = ff1.compute(gpos=gposnn1)
            gposnn2 = np.zeros(ff2.system.pos.shape, float)
            ff2.nlist.neighs[0] = (0, 1, d, 0.0, 0.0, d, 0, 0, 0)
            energy2 = ff2.compute(gpos=gposnn2)
            toplot.append([d,energy0,energy1,energy2])
#TODO
# I use the original line of code of Jelle
            row = [d, energy0-energy1+energy2, gposnn0[0,2]-gposnn1[0,2]+gposnn2[0,2]]
            #row = [d, energy2, gposnn2[0,2]]
            energies.append( row )
        energies = np.asarray(energies)
        toplot = np.asarray(toplot)
#TODO
# I outcommented every call to matplotlib        
        #plt.clf()
        #plt.plot(toplot[:,0]/angstrom,toplot[:,1]/kjmol,label='ei')
        #plt.plot(toplot[:,0]/angstrom,toplot[:,2]/kjmol,label='ei_noradii')
        #plt.plot(toplot[:,0]/angstrom,toplot[:,3]/kjmol,label='mm3')
        #plt.legend(loc='best')
        #plt.ylim([-1000.0,1000.0])

        ffa0 = ff.system.ffatypes[ff.system.ffatype_ids[number[0]]]
        ffa1 = ff.system.ffatypes[ff.system.ffatype_ids[number[1]]]
        if ffa0>ffa1:
            name = '%s-%s' % (ffa0,ffa1)
        else:
            name = '%s-%s' % (ffa1,ffa0)
        #plt.savefig('plots/%s.png'%name)
        ftab.write("%s\nN %d R %f %f\n\n" % (name, energies.shape[0], distances[0], distances[-1]))
        for irow, row in enumerate(energies):
            ftab.write("%05d %+13.8f %+21.12f %+21.12f\n" % (irow+1, row[0], row[1], row[2]))
        print name
        #break

if __name__=='__main__':
    chk_list=glob.glob('*.chk')
    if len(chk_list)==1:
        system = System.from_file(chk_list[0])
    else:
        raise Exception("Zero or more than one chk files, aborting (I don't know which one to pick)")
        
    log.set_level(log.silent)
    ff = ForceField.generate(system,'ff-lammps/pars.txt',rcut=15*angstrom, alpha_scale=3.2, gcut_scale=1.5, smooth_ei=True)
    log.set_level(log.low)
    write_lammps_table_jelle(ff)
