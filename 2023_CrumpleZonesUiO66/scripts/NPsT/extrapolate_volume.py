#! /usr/bin/env python
#Sven Rogge, Ruben Demuynck

from molmod.units import *
from yaff import *
import numpy as np
import h5py
import os

log.set_level(0)

## --------------- INPUT SECTION ----------------------- ##

# specify number of unit cells in the h5 file
no_uc = 8
# specify the chk to read the initial structure from
original_volume = 7850
fname = './init_%s.chk' %original_volume
# specify start, end and step of the volume grid (in Angstrom^3 per unit cell)
vol_start = 7000
vol_stop = 7825
vol_step = 25
# specify the directory to store the structure files
direc1='./structures_interpolated_from_%s/' %original_volume


## ----------------- MAIN ---------------------------- ##

# make directory to store the structure files
if not os.path.exists(direc1): os.makedirs(direc1)

# read data from original chk file
sys_old = System.from_file(fname)

# extract relevant data
numbers = sys_old.numbers
ffatypes = sys_old.ffatypes
ffatype_ids = sys_old.ffatype_ids
bonds = sys_old.bonds
masses = sys_old.masses
charges = sys_old.charges
pos_old = sys_old.pos
rvecs_old = sys_old.cell.rvecs
vol_old = np.linalg.det(rvecs_old)


for j in np.arange(vol_start, vol_stop+vol_step, vol_step):
    # determine scaling factor:
    volume = no_uc*j*angstrom**3
    scale = np.power(volume/vol_old, 1.0/3.0)
    print('rescaling with factor %.3f' %scale)
    pos = pos_old * scale
    rvecs = rvecs_old * scale
    sys = System(numbers, pos, ffatypes=ffatypes, ffatype_ids=ffatype_ids, bonds=bonds, rvecs=rvecs,
                 charges=charges, masses=masses)
    sys.to_file(str(direc1)+'init_%.0f.chk' % (volume / no_uc / angstrom ** 3))
