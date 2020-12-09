#!/usr/bin/env python

import numpy as np
from yaff import System
from molmod.units import angstrom

numbers = np.array([8,1,1,99],dtype=int)
pos = np.array([[ 0.0   , 0.0,     0.0],
                [ 0.7570, 0.5859,  0.0],
                [ -0.7570, 0.5859,  0.0],
                [ 0.0   , 0.15, 0.0],])*angstrom
ffatypes=['TO','TH','TH','TM']
bonds = [[1,0],[2,0],[3,0]]
system = System(numbers,pos,ffatypes=ffatypes,bonds=bonds)
system.to_file('water.chk')


numbers = np.array([7,7,99], dtype=int)
pos = np.zeros((3,3))
pos[0,2] = 0.55*angstrom
pos[2,2] = -0.55*angstrom
ffatypes = ['N2n','N2c','N2n']
bonds = [[1,0],[2,1]]
system = System(numbers,pos,ffatypes=ffatypes,bonds=bonds)
system.to_file('N2.chk')


numbers = np.array([18], dtype=int)
pos = np.zeros((1,3))
ffatypes = ['Ar']
bonds = []
system = System(numbers,pos,ffatypes=ffatypes,bonds=bonds)
system.to_file('argon.chk')

