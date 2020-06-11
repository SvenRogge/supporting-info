#! /usr/bin/env python

from molmod.io import dump_chk, load_chk
from molmod.molecules import Molecule
from molmod.molecular_graphs import MolecularGraph
from molmod.molecular_graphs import HasAtomNumber, HasNumNeighbors, HasNeighborNumbers, HasNeighbors
from molmod.graphs import CritAnd, CritNot
from molmod.periodic import periodic as pt
from molmod.units import angstrom

from yaff import System as YaffSystem
from yaff.log import log

import sys, os, numpy as np, cPickle

fn_chk = 'sample_temp.chk'
fn_fin = 'sample.chk'


# PERIODIC TABLE
''' RADII MOF-FF '''
radii = [(1,0.724),(6,1.163),(7,1.125),(8,1.118),(40,2.367)]

''' brick

ffatypes =  ['ZR', 'O_OX', 'O_OH', 'H_OH', 'O_CA', 'C_CA', 'C_PH', 'H_PH', 'C_PC']

afilters = [
    ('ZR' , HasAtomNumber(40)),
    ('O_OX', CritAnd(HasAtomNumber(8), HasNeighborNumbers(40, 40, 40))),
    ('O_OH',  CritAnd(HasAtomNumber(8), HasNeighborNumbers(40, 40, 40, 1))),
    ('H_OH',  CritAnd(HasAtomNumber(1), HasNeighborNumbers(8))),
    ('O_CA',  CritAnd(HasAtomNumber(8), HasNeighborNumbers(40, 6))),
    ('C_CA',  CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))),
    ('C_PH',  CritAnd(HasAtomNumber(6), HasNeighborNumbers(6, 6, 1))),
    ('H_PH',  CritAnd(HasAtomNumber(1), HasNeighborNumbers(6))),
    ('C_PC',  CritAnd(HasAtomNumber(6), HasNeighborNumbers(6, 6, 6))),
  ]

'''

''' linker 2 

ffatypes =  ['ZR', 'O_OX', 'O_OH', 'H_OH', 'O_CA', 'C_CA', 'C_PH', 'H_PH', 'C_PC']

afilters = [
    ('ZR' , HasAtomNumber(40)),
    ('O_OX', CritAnd(HasAtomNumber(8), HasNeighborNumbers(40, 40, 40))),
    ('O_OH',  CritAnd(HasAtomNumber(8), HasNeighborNumbers(40, 40, 40, 1))),
    ('H_OH',  CritAnd(HasAtomNumber(1), HasNeighborNumbers(8))),
    ('O_CA',  CritAnd(HasAtomNumber(8), HasNeighborNumbers(40, 6))),
    ('C_CA',  CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))),
    ('C_PH',  CritAnd(HasAtomNumber(6), HasNeighborNumbers(6, 6, 1))),
    ('H_PH',  CritAnd(HasAtomNumber(1), HasNeighbors(CritAnd(HasAtomNumber(6), HasNeighborNumbers(6, 6, 1))))),
    ('C_PC',  CritAnd(HasAtomNumber(6), HasNeighborNumbers(6, 6, 6))),
    ('c_ca', CritAnd(HasAtomNumber(6), HasNeighborNumbers(8,8,1))),
    ('h_ca', CritAnd(HasAtomNumber(1), HasNeighbors(CritAnd(HasAtomNumber(6), HasNeighborNumbers(8,8,1))))),
  ]

'''

''' linker F2 

ffatypes =  ['ZR', 'O_OX', 'O_OH', 'H_OH', 'O_CA', 'C_CA', 'C_PH', 'H_PH', 'C_PC']

afilters = [
    ('ZR' , HasAtomNumber(40)),
    ('O_OX', CritAnd(HasAtomNumber(8), HasNeighborNumbers(40, 40, 40))),
    ('O_OH',  CritAnd(HasAtomNumber(8), HasNeighborNumbers(40, 40, 40, 1))),
    ('H_OH',  CritAnd(HasAtomNumber(1), HasNeighbors(HasNeighborNumbers(40,40,40,1)))),
    ('O_CA',  CritAnd(HasAtomNumber(8), HasNeighborNumbers(40, 6))),
    ('O_C', CritAnd(HasAtomNumber(8), HasNeighborNumbers(6))),
    ('O_CH', CritAnd(HasAtomNumber(8), HasNeighborNumbers(6,1))),
    ('C_O', CritAnd(HasAtomNumber(6), HasNeighbors(CritAnd(HasAtomNumber(8), HasNeighborNumbers(6)), CritAnd(HasAtomNumber(8), HasNeighborNumbers(6,1)), HasAtomNumber(6)))),
    ('C_CA',  CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))),
    ('C_C', CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighbors(CritAnd(HasAtomNumber(8), HasNeighborNumbers(6)), CritAnd(HasAtomNumber(8), HasNeighborNumbers(6,1)), HasAtomNumber(6)))))),
    ('H_O', CritAnd(HasAtomNumber(1), HasNeighbors(CritAnd(HasAtomNumber(8), HasNeighborNumbers(6,1))))),
    ('C_H', CritAnd(HasAtomNumber(6), HasNeighborNumbers(6,6,1))),
    ('H', CritAnd(HasAtomNumber(1), HasNeighbors(CritAnd(HasAtomNumber(6), HasNeighborNumbers(6,6,1))))),
    ('C_PC',  CritAnd(HasAtomNumber(6), HasNeighborNumbers(6, 6, 6))),
    ('cl_ca', CritAnd(HasAtomNumber(6), HasNeighborNumbers(8,8,1))),
    ('hl_ca', CritAnd(HasAtomNumber(1), HasNeighbors(CritAnd(HasAtomNumber(6), HasNeighborNumbers(8,8,1))))),
  ]
'''


''' linker 13 '''
ffatypes =  ['ZR', 'O_OX', 'O_OH', 'H_OH', 'O_CA', 'C_CA', 'C_AL', 'c_ca', 'h_ca']

afilters = [
    ('ZR' , HasAtomNumber(40)),
    ('O_OX', CritAnd(HasAtomNumber(8), HasNeighborNumbers(40, 40, 40))),
    ('O_OH',  CritAnd(HasAtomNumber(8), HasNeighborNumbers(40, 40, 40, 1))),
    ('H_OH',  CritAnd(HasAtomNumber(1), HasNeighborNumbers(8))),
    ('O_CA',  CritAnd(HasAtomNumber(8), HasNeighborNumbers(40, 6))),
    ('C_CA',  CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))),
    ('C1', CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(7), HasAtomNumber(7), CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))))),
    ('N1', CritAnd(HasAtomNumber(7), HasNeighbors(HasAtomNumber(7), CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(7), HasAtomNumber(7), CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))))))),
    ('N2', CritAnd(HasAtomNumber(7), HasNeighbors(HasAtomNumber(6), CritAnd(HasAtomNumber(7), HasNeighbors(HasAtomNumber(7), CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(7), HasAtomNumber(7), CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))))))))),
    ('C2', CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(7), CritAnd(HasAtomNumber(7), HasNeighbors(HasAtomNumber(6), CritAnd(HasAtomNumber(7), HasNeighbors(HasAtomNumber(7), CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(7), HasAtomNumber(7), CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))))))))))),
    ('C3', CritAnd(HasNeighbors(HasAtomNumber(7), HasAtomNumber(7), CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(7), CritAnd(HasAtomNumber(7), HasNeighbors(HasAtomNumber(6), CritAnd(HasAtomNumber(7), HasNeighbors(HasAtomNumber(7), CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(7), HasAtomNumber(7), CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))))))))))))),
    ('N4', CritAnd(HasNeighbors(HasAtomNumber(7), CritAnd(HasNeighbors(HasAtomNumber(7), HasAtomNumber(7), CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(7), CritAnd(HasAtomNumber(7), HasNeighbors(HasAtomNumber(6), CritAnd(HasAtomNumber(7), HasNeighbors(HasAtomNumber(7), CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(7), HasAtomNumber(7), CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))))))))))))))),
    ('cl_ca', CritAnd(HasAtomNumber(6), HasNeighborNumbers(8,8,1))),
    ('hl_ca', CritAnd(HasAtomNumber(1), HasNeighbors(CritAnd(HasAtomNumber(6), HasNeighborNumbers(8,8,1))))),
  ]
  
''' linker 7 '''
ffatypes =  ['ZR', 'O_OX', 'O_OH', 'H_OH', 'O_CA', 'C_CA', 'C_AL', 'c_ca', 'h_ca']

afilters = [
    ('ZR' , HasAtomNumber(40)),
    ('O_OX', CritAnd(HasAtomNumber(8), HasNeighborNumbers(40, 40, 40))),
    ('O_OH',  CritAnd(HasAtomNumber(8), HasNeighborNumbers(40, 40, 40, 1))),
    ('H_OH',  CritAnd(HasAtomNumber(1), HasNeighborNumbers(8))),
    ('O_CA',  CritAnd(HasAtomNumber(8), HasNeighborNumbers(40, 6))),
    ('C_CA',  CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))),
    ('C1', CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(7), HasAtomNumber(7), CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))))),
    ('N1', CritAnd(HasAtomNumber(7), HasNeighbors(HasAtomNumber(7), CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(7), HasAtomNumber(7), CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))))))),
    ('N2', CritAnd(HasAtomNumber(7), HasNeighbors(HasAtomNumber(6), CritAnd(HasAtomNumber(7), HasNeighbors(HasAtomNumber(7), CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(7), HasAtomNumber(7), CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))))))))),
    ('C2', CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(7), CritAnd(HasAtomNumber(7), HasNeighbors(HasAtomNumber(6), CritAnd(HasAtomNumber(7), HasNeighbors(HasAtomNumber(7), CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(7), HasAtomNumber(7), CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))))))))))),    
    ('cl_ca', CritAnd(HasAtomNumber(6), HasNeighborNumbers(8,8,1))),
    ('hl_ca', CritAnd(HasAtomNumber(1), HasNeighbors(CritAnd(HasAtomNumber(6), HasNeighborNumbers(8,8,1))))),
  ]
  
''' linker 3 '''
ffatypes =  ['ZR', 'O_OX', 'O_OH', 'H_OH', 'O_CA', 'C_CA', 'C_AL', 'c_ca', 'h_ca']

afilters = [
    ('ZR' , HasAtomNumber(40)),
    ('O_OX', CritAnd(HasAtomNumber(8), HasNeighborNumbers(40, 40, 40))),
    ('O_OH',  CritAnd(HasAtomNumber(8), HasNeighborNumbers(40, 40, 40, 1))),
    ('H_OH',  CritAnd(HasAtomNumber(1), HasNeighborNumbers(8))),
    ('O_CA',  CritAnd(HasAtomNumber(8), HasNeighborNumbers(40, 6))),
    ('C_CA',  CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))),
    ('C1', CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(7), HasAtomNumber(7), CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))))),
    ('N1', CritAnd(HasAtomNumber(7), HasNeighbors(HasAtomNumber(7), CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(7), HasAtomNumber(7), CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))))))),
    ('cl_ca', CritAnd(HasAtomNumber(6), HasNeighborNumbers(8,8,1))),
    ('hl_ca', CritAnd(HasAtomNumber(1), HasNeighbors(CritAnd(HasAtomNumber(6), HasNeighborNumbers(8,8,1))))),
  ]

''' linker 1 

ffatypes =  ['ZR', 'O_OX', 'O_OH', 'H_OH', 'O_CA', 'C_CA', 'C_AL', 'c_ca', 'h_ca']

afilters = [
    ('ZR' , HasAtomNumber(40)),
    ('O_OX', CritAnd(HasAtomNumber(8), HasNeighborNumbers(40, 40, 40))),
    ('O_OH',  CritAnd(HasAtomNumber(8), HasNeighborNumbers(40, 40, 40, 1))),
    ('H_OH',  CritAnd(HasAtomNumber(1), HasNeighborNumbers(8))),
    ('O_CA',  CritAnd(HasAtomNumber(8), HasNeighborNumbers(40, 6))),
    ('C_CA',  CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))),
    ('C_AL',  CritAnd(HasAtomNumber(6), HasNeighborNumbers(6, 6))),
    ('c_ca', CritAnd(HasAtomNumber(6), HasNeighborNumbers(8,8,1))),
    ('h_ca', CritAnd(HasAtomNumber(1), HasNeighbors(CritAnd(HasAtomNumber(6), HasNeighborNumbers(8,8,1))))),
  ]
  
'''
  
''' linker 4 
  
ffatypes =  ['ZR', 'O_OX', 'O_OH', 'H_OH', 'O_CA', 'C_CA', 'C_AL', 'c_ca', 'h_ca']

afilters = [
    ('ZR' , HasAtomNumber(40)),
    ('O_OX', CritAnd(HasAtomNumber(8), HasNeighborNumbers(40, 40, 40))),
    ('O_OH',  CritAnd(HasAtomNumber(8), HasNeighborNumbers(40, 40, 40, 1))),
    ('H_OH',  CritAnd(HasAtomNumber(1), HasNeighborNumbers(8))),
    ('O_CA',  CritAnd(HasAtomNumber(8), HasNeighborNumbers(40, 6))),
    ('C_CA',  CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))),
    ('C_PC',  CritAnd(HasAtomNumber(6), HasNeighborNumbers(6, 6), HasNeighbors(HasAtomNumber(6),HasNeighborNumbers(8, 8, 6)))),
    ('C_AL',  CritAnd(HasAtomNumber(6), HasNeighborNumbers(6, 6), HasNeighbors(HasAtomNumber(6),HasNeighborNumbers(6, 6)))),
    ('c_ca', CritAnd(HasAtomNumber(6), HasNeighborNumbers(8,8,1))),
    ('h_ca', CritAnd(HasAtomNumber(1), HasNeighbors(CritAnd(HasAtomNumber(6), HasNeighborNumbers(8,8,1))))),
  ]
  
'''

''' linker 5 

ffatypes =  ['ZR', 'O_OX', 'O_OH', 'H_OH', 'O_CA', 'C_CA', 'C1', 'C6', 'C5', 'C2', 'C3', 'C4', 'c_ca', 'h_ca']

afilters = [
    ('ZR' , HasAtomNumber(40)),
    ('O_OX', CritAnd(HasAtomNumber(8), HasNeighborNumbers(40, 40, 40))),
    ('O_OH',  CritAnd(HasAtomNumber(8), HasNeighborNumbers(40, 40, 40, 1))),
    ('H_OH',  CritAnd(HasAtomNumber(1), HasNeighborNumbers(8))),
    ('O_CA',  CritAnd(HasAtomNumber(8), HasNeighborNumbers(40, 6))),
    ('C_CA',  CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))),
    ('C1',  CritAnd(HasAtomNumber(6), HasNeighborNumbers(6, 6, 6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))))),
    ('C6', CritAnd(HasAtomNumber(6), HasNeighborNumbers(6, 6), HasNeighbors(HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))))),
    ('C5', CritAnd(HasAtomNumber(6), HasNeighborNumbers(6,6), HasNeighbors(CritAnd(HasAtomNumber(6), HasNeighborNumbers(6, 6), HasNeighbors(HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6)))), HasAtomNumber(6)))),
    ('C2', CritAnd(HasAtomNumber(6), HasNeighborNumbers(6,6,1), HasNeighbors(HasAtomNumber(1), HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighborNumbers(6, 6, 6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))))))),
    ('C3', CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(1), CritAnd(HasAtomNumber(6), HasNeighborNumbers(6,6,1), HasNeighbors(HasAtomNumber(1), HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighborNumbers(6, 6, 6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))))))))),
    ('C4', CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(1), HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighborNumbers(6,6,1), HasNeighbors(HasAtomNumber(1), HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighborNumbers(6, 6, 6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))))))))))),   
    ('H2', CritAnd(HasAtomNumber(1), HasNeighbors(CritAnd(HasAtomNumber(6), HasNeighborNumbers(6,6,1), HasNeighbors(HasAtomNumber(1), HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighborNumbers(6, 6, 6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))))))))),
    ('H3', CritAnd(HasAtomNumber(1), HasNeighbors(CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(1), CritAnd(HasAtomNumber(6), HasNeighborNumbers(6,6,1), HasNeighbors(HasAtomNumber(1), HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighborNumbers(6, 6, 6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))))))))))), 
    ('c_ca', CritAnd(HasAtomNumber(6), HasNeighborNumbers(8,8,1))),
    ('h_ca', CritAnd(HasAtomNumber(1), HasNeighbors(CritAnd(HasAtomNumber(6), HasNeighborNumbers(8,8,1))))),
  ]

'''


''' linker 8 

ffatypes =  ['ZR', 'O_OX', 'O_OH', 'H_OH', 'O_CA', 'C_CA', 'C1', 'C6', 'C5', 'C2', 'C3', 'C4', 'c_ca', 'h_ca']

afilters = [
    ('ZR' , HasAtomNumber(40)),
    ('O_OX', CritAnd(HasAtomNumber(8), HasNeighborNumbers(40, 40, 40))),
    ('O_OH',  CritAnd(HasAtomNumber(8), HasNeighborNumbers(40, 40, 40, 1))),
    ('H_OH',  CritAnd(HasAtomNumber(1), HasNeighborNumbers(8))),
    ('O_CA',  CritAnd(HasAtomNumber(8), HasNeighborNumbers(40, 6))),
    ('C_CA',  CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))),
    ('C1', CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))))),
    ('C2', CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))))))),
    ('C3', CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))))))))),
    ('c_ca', CritAnd(HasAtomNumber(6), HasNeighborNumbers(8,8,1))),
    ('h_ca', CritAnd(HasAtomNumber(1), HasNeighbors(CritAnd(HasAtomNumber(6), HasNeighborNumbers(8,8,1))))),
    ]    
'''   

''' linker 9 

ffatypes =  ['ZR', 'O_OX', 'O_OH', 'H_OH', 'O_CA', 'C_CA', 'C_AL', 'c_ca', 'h_ca']

afilters = [
    ('ZR' , HasAtomNumber(40)),
    ('O_OX', CritAnd(HasAtomNumber(8), HasNeighborNumbers(40, 40, 40))),
    ('O_OH',  CritAnd(HasAtomNumber(8), HasNeighborNumbers(40, 40, 40, 1))),
    ('H_OH',  CritAnd(HasAtomNumber(1), HasNeighborNumbers(8))),
    ('O_CA',  CritAnd(HasAtomNumber(8), HasNeighborNumbers(40, 6))),
    ('C_CA',  CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))),
    ('C1', CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(7), HasAtomNumber(7), CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))))),
    ('N1', CritAnd(HasAtomNumber(7), HasNeighbors(HasAtomNumber(7), CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(7), HasAtomNumber(7), CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))))))),
    ('N2', CritAnd(HasAtomNumber(7), HasNeighbors(HasAtomNumber(6), CritAnd(HasAtomNumber(7), HasNeighbors(HasAtomNumber(7), CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(7), HasAtomNumber(7), CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))))))))),
    ('C2', CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(7), CritAnd(HasAtomNumber(7), HasNeighbors(HasAtomNumber(6), CritAnd(HasAtomNumber(7), HasNeighbors(HasAtomNumber(7), CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(7), HasAtomNumber(7), CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))))))))))),
    ('C3', CritAnd(HasNeighbors(HasAtomNumber(7), HasAtomNumber(7), CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(7), CritAnd(HasAtomNumber(7), HasNeighbors(HasAtomNumber(6), CritAnd(HasAtomNumber(7), HasNeighbors(HasAtomNumber(7), CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(7), HasAtomNumber(7), CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))))))))))))),
    ('N4', CritAnd(HasNeighbors(HasAtomNumber(7), CritAnd(HasNeighbors(HasAtomNumber(7), HasAtomNumber(7), CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(7), CritAnd(HasAtomNumber(7), HasNeighbors(HasAtomNumber(6), CritAnd(HasAtomNumber(7), HasNeighbors(HasAtomNumber(7), CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(7), HasAtomNumber(7), CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6)))))))))))))))
    ] HasAtomNumber(8), HasNeighborNumbers(40, 6))),
    ('C_CA',  CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))),
    ('C1', CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))))),
    ('C2', CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))))))),
    ('C3', CritAnd(HasAtomNumber(6), HasNeighborNumbers(6,6,6))),
    ('C4', CritAnd(HasAtomNumber(6), HasNeighborNumbers(6,6,1))),
    ('H4', CritAnd(HasAtomNumber(1), HasNeighbors(CritAnd(HasAtomNumber(6), HasNeighborNumbers(6,6,1))))),
    ('c_ca', CritAnd(HasAtomNumber(6), HasNeighborNumbers(8,8,1))),
    ('h_ca', CritAnd(HasAtomNumber(1), HasNeighbors(CritAnd(HasAtomNumber(6), HasNeighborNumbers(8,8,1))))),
  ]  
    
'''

''' linker 10 
ffatypes =  ['ZR', 'O_OX', 'O_OH', 'H_OH', 'O_CA', 'C_CA', 'C_AL', 'c_ca', 'h_ca']

afilters = [
    ('ZR' , HasAtomNumber(40)),
    ('O_OX', CritAnd(HasAtomNumber(8), HasNeighborNumbers(40, 40, 40))),
    ('O_OH',  CritAnd(HasAtomNumber(8), HasNeighborNumbers(40, 40, 40, 1))),
    ('H_OH',  CritAnd(HasAtomNumber(1), HasNeighborNumbers(8))),
    ('O_CA',  CritAnd(HasAtomNumber(8), HasNeighborNumbers(40, 6))),
    ('C_CA',  CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))),
    ('C8', CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))))),
    ('C7', CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))))))),
    ('C6', CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))))))))),
    ('C5', CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))))))))))),
    ('C1', CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))))),
    ('C2', CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(1), CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))))))),
    ('C3', CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(1), CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(1), CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))))))))),
    ('C4', CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(1), CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(1), CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))))))))))),
    ('H2', CritAnd(HasAtomNumber(1), HasNeighbors(CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(1), CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))))))))),
    ('H3', CritAnd(HasAtomNumber(1), HasNeighbors(CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(1), CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(1), CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))))))))))),
    ('c_ca', CritAnd(HasAtomNumber(6), HasNeighborNumbers(8,8,1))),
    ('h_ca', CritAnd(HasAtomNumber(1), HasNeighbors(CritAnd(HasAtomNumber(6), HasNeighborNumbers(8,8,1))))),   
    ]
    
'''
    
    

''' linker 11 
ffatypes =  ['ZR', 'O_OX', 'O_OH', 'H_OH', 'O_CA', 'C_CA', 'C_AL', 'c_ca', 'h_ca']

afilters = [
    ('ZR' , HasAtomNumber(40)),
    ('O_OX', CritAnd(HasAtomNumber(8), HasNeighborNumbers(40, 40, 40))),
    ('O_OH',  CritAnd(HasAtomNumber(8), HasNeighborNumbers(40, 40, 40, 1))),
    ('H_OH',  CritAnd(HasAtomNumber(1), HasNeighborNumbers(8))),
    ('O_CA',  CritAnd(HasAtomNumber(8), HasNeighborNumbers(40, 6))),
    ('C_CA',  CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))),
    ('C5',  CritAnd(HasAtomNumber(6), HasNeighborNumbers(6, 6))),
    ('C4', CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighborNumbers(6, 6))))),
    ('C1', CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))))),
    ('C2', CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(1),HasAtomNumber(6),CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))))))),
    ('C3', CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(1),HasAtomNumber(6),CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(1),HasAtomNumber(6),CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))))))))),
    ('H2', CritAnd(HasAtomNumber(1), HasNeighbors(CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(1),HasAtomNumber(6),CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))))))))),
    ('H3', CritAnd(HasAtomNumber(1), HasNeighbors(CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(1),HasAtomNumber(6),CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(1),HasAtomNumber(6),CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))))))))))),
    ('c_ca', CritAnd(HasAtomNumber(6), HasNeighborNumbers(8,8,1))),
    ('h_ca', CritAnd(HasAtomNumber(1), HasNeighbors(CritAnd(HasAtomNumber(6), HasNeighborNumbers(8,8,1))))),
  ]
'''

''' linker 14 
ffatypes =  ['ZR', 'O_OX', 'O_OH', 'H_OH', 'O_CA', 'C_CA', 'C_AL', 'c_ca', 'h_ca']

afilters = [
    ('ZR' , HasAtomNumber(40)),
    ('O_OX', CritAnd(HasAtomNumber(8), HasNeighborNumbers(40, 40, 40))),
    ('O_OH',  CritAnd(HasAtomNumber(8), HasNeighborNumbers(40, 40, 40, 1))),
    ('H_OH',  CritAnd(HasAtomNumber(1), HasNeighborNumbers(8))),
    ('O_CA',  CritAnd(HasAtomNumber(8), HasNeighborNumbers(40, 6))),
    ('C_CA',  CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))),
    ('C1', CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))))),
    ('C2', CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(1), CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))))))),
    ('C3', CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(1), CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(1), CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))))))))),
    ('C4', CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(1), CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(1), CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))))))))))),
    ('C5', CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(1), CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(1), CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))))))))))))),
    ('C6', CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(1), CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(1), CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))))))))))))))),
    ('H2', CritAnd(HasAtomNumber(1), HasNeighbors(CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(1), CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))))))))),
    ('H3', CritAnd(HasAtomNumber(1), HasNeighbors(CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(1), CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(1), CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))))))))))),
    ('c_ca', CritAnd(HasAtomNumber(6), HasNeighborNumbers(8,8,1))),
    ('h_ca', CritAnd(HasAtomNumber(1), HasNeighbors(CritAnd(HasAtomNumber(6), HasNeighborNumbers(8,8,1))))),
    ]    
 
'''



'''UiO-67
ffatypes = ['ZR', 'O_OX', 'O_OH', 'H_OH', 'O_CA', 'C_CA', 'C1', 'C2', 'C3', 'C4', 'H2', 'H3', 'c_ca', 'h_ca']

afilters = [
    ('ZR' , HasAtomNumber(40)),
    ('O_OX', CritAnd(HasAtomNumber(8), HasNeighborNumbers(40, 40, 40))),
    ('O_OH',  CritAnd(HasAtomNumber(8), HasNeighborNumbers(40, 40, 40, 1))),
    ('H_OH',  CritAnd(HasAtomNumber(1), HasNeighborNumbers(8))),
    ('O_CA',  CritAnd(HasAtomNumber(8), HasNeighborNumbers(40, 6))),
    ('C_CA',  CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))),
    ('C1',  CritAnd(HasAtomNumber(6), HasNeighbors(HasNeighborNumbers(6, 6, 1), HasNeighborNumbers(6, 6, 1), HasNeighborNumbers(6, 8, 8)))),
    ('C2', CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(1), HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighbors(HasNeighborNumbers(6, 6, 1), HasNeighborNumbers(6, 6, 1), HasNeighborNumbers(6, 8, 8)))))),
    ('C3', CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(1), HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighbors(HasNeighborNumbers(6,6,6), HasNeighborNumbers(6,6,1), HasNeighborNumbers(6,6,1)))))),
    ('C4', CritAnd(HasAtomNumber(6), HasNeighbors(HasNeighborNumbers(6,6,6), HasNeighborNumbers(6,6,1), HasNeighborNumbers(6,6,1)))),
    ('H2', CritAnd(HasAtomNumber(1), HasNeighbors(CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(1), HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighbors(HasNeighborNumbers(6, 6, 1), HasNeighborNumbers(6, 6, 1), HasNeighborNumbers(6, 8, 8)))))))),
    ('H3', CritAnd(HasAtomNumber(1), HasNeighbors(CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(1), HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighbors(HasNeighborNumbers(6,6,6), HasNeighborNumbers(6,6,1), HasNeighborNumbers(6,6,1)))))))),
    ('c_ca', CritAnd(HasAtomNumber(6), HasNeighborNumbers(8,8,1))),
    ('h_ca', CritAnd(HasAtomNumber(1), HasNeighbors(CritAnd(HasAtomNumber(6), HasNeighborNumbers(8,8,1)))))
  ]
'''

'''UiO-68
ffatypes = ['ZR', 'O_OX', 'O_OH', 'H_OH', 'O_CA', 'C_CA', 'C1', 'C2', 'C3', 'C4', 'H1', 'H2']

afilters = [
    ('ZR' , HasAtomNumber(40)),
    ('O_OX', CritAnd(HasAtomNumber(8), HasNeighborNumbers(40, 40, 40))),
    ('O_OH',  CritAnd(HasAtomNumber(8), HasNeighborNumbers(40, 40, 40, 1))),
    ('H_OH',  CritAnd(HasAtomNumber(1), HasNeighborNumbers(8))),
    ('O_CA',  CritAnd(HasAtomNumber(8), HasNeighborNumbers(40, 6))),
    ('C_CA',  CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))),
    ('C1',  CritAnd(HasAtomNumber(6), HasNeighbors(HasNeighborNumbers(6, 6, 1), HasNeighborNumbers(6, 6, 1), HasNeighborNumbers(6, 8, 8)))),
     ('C1', CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))))),
    ('C2', CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(1),HasAtomNumber(6),CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))))))),
    ('C3', CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(1),HasAtomNumber(6),CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(1),HasAtomNumber(6),CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))))))))),
    ('C4', CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(1),HasAtomNumber(6),CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(1),HasAtomNumber(6),CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))))))))))),
    ('C5', CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(1),HasAtomNumber(6),CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(1),HasAtomNumber(6),CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))))))))))))),
    ('C6', CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(1), HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(1),HasAtomNumber(6),CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(1),HasAtomNumber(6),CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))))))))))))))),
    ('H2', CritAnd(HasAtomNumber(1), HasNeighbors(CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(1),HasAtomNumber(6),CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))))))))),
    ('H3', CritAnd(HasAtomNumber(1), HasNeighbors(CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(1),HasAtomNumber(6),CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(1),HasAtomNumber(6),CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))))))))))),
    ('H6', CritAnd(HasAtomNumber(1), HasNeighbors(CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(1), HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(1),HasAtomNumber(6),CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(1),HasAtomNumber(6),CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(6), HasAtomNumber(6), CritAnd(HasAtomNumber(6), HasNeighborNumbers(8, 8, 6))))))))))))))))),    
    ('c_ca', CritAnd(HasAtomNumber(6), HasNeighborNumbers(8,8,1))),
    ('h_ca', CritAnd(HasAtomNumber(1), HasNeighbors(CritAnd(HasAtomNumber(6), HasNeighborNumbers(8,8,1)))))  
  ]
  
'''


def get_atypes(yaffsystem,ffatypes):
    graph = MolecularGraph(yaffsystem.bonds, yaffsystem.numbers)
    ffatype_ids = [0]*len(yaffsystem.numbers)
    ffatypes_2 = [0]*len(yaffsystem.numbers)
    test = [False]*len(yaffsystem.numbers)
    teller = -1
    for ffatype, filter in afilters:
        print ffatype
        teller = teller + 1 
        for iatom, number in enumerate(yaffsystem.numbers):
            if filter(iatom, graph) and not test[iatom]:
                ffatype_ids[iatom] = teller
                test[iatom] = True
                ffatypes_2[iatom] = ffatype
                print ffatype, iatom
    for iatom, number in enumerate(yaffsystem.numbers):
        if not test[iatom]:
            print 'No atom type found for atom %i(%s)' %(iatom, pt[number].symbol)
            print 'This atom has neighbors:'
            for neighbor in graph.neighbors[iatom]:
                print '  %i (%s)' %(neighbor, pt[yaffsystem.numbers[neighbor]].symbol)
    return ffatype_ids, ffatypes_2


'''
def get_atypes(yaffsystem,ffatypes):
    ffatype_ids = [0]*len(yaffsystem.numbers)
    unique_atypes = []
    unique_atypes2 = []
    teller = 0
    for atype in ffatypes:                        
        if atype not in unique_atypes:
            unique_atypes.append(atype)
            unique_atypes2.append(teller)
            teller += 1
    for I, atom in enumerate(ffatypes):
        for J, fftype in enumerate(unique_atypes):
            if fftype is atom:
                ffatype_ids[I] = unique_atypes2[J]
    return ffatype_ids, ffatypes
'''

def read_sample(fn,fn_fin,ffatypes):
    #log.set_level(log.silent)
    yaffsystem = YaffSystem.from_file(fn)    
    
    # Manually adjust if necessary!!!                                                                          !    !   !   !   !   !
    GaGabond = {(40,40):1.0*angstrom,(40,6):0.5*angstrom,(1,40):0.5*angstrom}
    #VVbond = {(22,22):0.01, (13,13):0.01, (24,24):0.01, (31,31):0.01, (21,21):0.01}
    #AlAlbond = {(13,13):0.1}
    
    yaffsystem.detect_bonds(GaGabond)
    
    #yaffsystem.detect_bonds()
    
    radii_list = []
    for I,number in enumerate(yaffsystem.numbers):
        test = False
        for J,element in enumerate(radii):
            if int(element[0]) == number:
                radii_list.append(float(element[1])*angstrom)
                test = True
        if not test:
            print 'NO RADIUS FOR ELEMENT %s'%element[0]
    
    # Set atomtypes
    ffatype_ids, ffatypes_2 = get_atypes(yaffsystem, ffatypes)
    # Finish
    temp = load_chk(fn)
    dump_chk(fn_fin, {
    'numbers': yaffsystem.numbers,
    'pos': yaffsystem.pos,
    'ffatypes': ffatypes_2,
    'bonds': yaffsystem.bonds,
    #'rvecs': yaffsystem.cell.rvecs,
    #'masses': yaffsystem.masses,
    #'hessian': temp['hessian'],
    #'gradient': temp['gradient'],
    'radii': radii_list
    })

def main():
    #READ SAMPLE
    read_sample(fn_chk,fn_fin,ffatypes)                   
    
if __name__=='__main__':
    main()
