#!/usr/bin/env python

from yaff import *
import h5py
from molmod.io.xyz import XYZReader
from molmod.units import *
import numpy as np
import sys
import matplotlib.pyplot as pt

rcut = 10*angstrom
rspacing = 0.02*angstrom
nimage = 2

def make_rdf(rcut, rspacing, f, group0, group1, nimage, cellpath, name):
    rdf = RDF(rcut=rcut, rspacing=rspacing, f=f, select0=group0, select1=group1, nimage=nimage, cellpath=cellpath)
    rdf.compute_iteration()
    rdf.compute_derived()
    g = open('RDF_'+name+'.dat', 'w')
    g.write('#Snap\tDistance\tRDF\tCRDF\n')
    for i in range(len(rdf.d)):
        g.write(str(rdf.d[i]/angstrom)+"\t"+str(rdf.rdf[i])+"\t"+str(rdf.rdf[:i+1].sum()*rspacing)+"\n")
    g.close()


with h5py.File('rdf.h5', 'a') as f:
    numbers = np.array(f['system/numbers'])
    cellpath = 'trajectory/cell'

    id_dummy = np.where(numbers==84)[0]
    id_o = np.where(numbers==8)[0]

    make_rdf(rcut, rspacing, f, id_dummy, id_o, nimage, cellpath, 'O')
    make_rdf(rcut, rspacing, f, id_dummy, None, nimage, cellpath, '6MR')

rdf_o = np.loadtxt('RDF_O.dat')
rdf_6MR = np.loadtxt('RDF_6MR.dat')

pt.clf()
fig, ax = pt.subplots()

max_rdf6MR = np.amax(rdf_6MR[:,1])
ax.axvline(rdf_6MR[0,0],color='#fec44f',linewidth=19.5*rspacing,alpha=1)
for i in range(len(rdf_6MR[:,0])):
    ax.axvline(rdf_6MR[i,0],color='#fec44f',linewidth=19.5*rspacing,alpha=rdf_6MR[i,1]/max_rdf6MR)

pt.plot(rdf_o[:,0], rdf_o[:,1], color='0.3')

pt.xlim([0,10])
pt.xlabel('Distance from center of the sixring [A]')
pt.ylabel('Radial distribution function [1/A]')
pt.savefig('RDF.pdf', format='pdf', bbox_inches = 'tight')
pt.savefig('RDF.svg', format='svg', bbox_inches = 'tight')
