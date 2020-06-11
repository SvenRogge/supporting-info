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
no_uc = 1
# specify the directory to store the structure files
direc1='../structures/'
# specify the h5 file to read the initial structure from
fname = './traj_1.h5'
# specify start, end and step of the volume grid (in Angstrom^3 per unit cell)
vol_start = 5500
vol_stop = 6500
vol_step = 25
# specify the allowed volume deviation (in Angstrom^3 per unit cell)
vol_eps = 12


## ----------------- MAIN ---------------------------- ##

def determine_index(globalList, volume, volume_eps, no_uc):
    '''
        Returns the index where the volume 'volume' can be found
        in the sorted array 'globalList'. If no volume is found so
        that |volume_found - volume| < volume_eps, -1 is returned.
    '''

    # find the correct index using binary search algorithm
    indexGL = np.searchsorted(globalList[:,0], volume)
    # verify which neighbour is the closest
    if indexGL <> 0 and indexGL <> globalList.shape[0]:
        volume_lower = globalList[indexGL-1,0]
        volume_larger = globalList[indexGL,0]
        if volume - volume_lower < volume_larger - volume: indexGL -= 1
    if indexGL == globalList.shape[0]: indexGL -= 1
    # verify whether the suggestion is within volume_eps
    if np.abs(globalList[indexGL,0] - volume) <= volume_eps: 
        index = globalList[indexGL,1]
        scale = np.power(volume/globalList[indexGL,0], 1.0/3.0)
        print 'I have found a matching structure'
        print 'You were asking for volume %.0f A^3 and I found one with volume %.0f A^3' % (volume/angstrom**3/no_uc, globalList[indexGL,0]/angstrom**3/no_uc)
        return index, scale
    else:
        print 'I did not find a matching structure for volume %.0f A^3' % (volume/angstrom**3/no_uc)
        print 'You trusted me, and I failed you...'
        return -1, 1

# make directory to store the structure files
if not os.path.exists(direc1): os.makedirs(direc1)

# read data from h5 file
f = h5py.File(fname)

# data which is constant for all volumes
numbers = np.array(f['system/numbers'])
ffatypes = f['system/ffatypes']
ffatype_ids = np.array(f['system/ffatype_ids'])
bonds = np.array(f['system/bonds'])
masses = np.array(f['system/masses'])
charges = np.array(f['system/charges'])

# initialize the data arrays
volumes = f['trajectory/volume']
globalList = np.zeros((len(volumes),2))
globalList[:,0] = volumes
globalList[:,1] = np.arange(0, len(volumes))

# sort the data according to the volume
ind = np.lexsort((globalList[:,1], globalList[:,0]))
globalList[:] = globalList[ind]

# volume deviation for the super cell
volume_eps = no_uc*vol_eps*angstrom**3

for j in np.arange(vol_start, vol_stop+vol_step, vol_step):
    # super cell volume to look for
    volume = no_uc*j*angstrom**3
    index, scale = determine_index(globalList, volume, volume_eps, no_uc)
    # if super cell volume is found: write initial structure file
    if index <> -1:
        pos = f['trajectory/pos'][index, :, :] * scale
        rvecs = f['trajectory/cell'][index, :, :] * scale
        sys = System(numbers, pos, ffatypes=ffatypes, ffatype_ids=ffatype_ids, bonds=bonds, rvecs=rvecs,
                     charges=charges, masses=masses)
        sys.to_file(str(direc1)+'init_%.0f.chk' % (volume / no_uc / angstrom ** 3))
