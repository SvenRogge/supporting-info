#!/usr/bin/env python

import numpy as np
from yaff import *
from molmod.units import angstrom
import itertools

chk_file = 'init.chk'
ff_file = 'pars.txt'

# -------PART 1: ADAPTING THE SYSTEM OBJECT----- #
system = System.from_file(chk_file)

ffatype_1 = "O_CA"
ffatype_2 = "O_OX"
ffatype_3 = "O_OH"

idx_ce = 72 # for hafnium; for cerium: 58
idx_zr = 40

ffidx_1 = -1
ffidx_2 = -1
ffidx_3 = -1
ffidx_max = -1

ffatype_ids_copy = system.ffatype_ids.copy()

for idx, ffatype in enumerate(system.ffatypes):
    if(ffatype == ffatype_1): ffidx_1 = idx       # determine index of O_CA
    elif(ffatype == ffatype_2): ffidx_2 = idx     # determine index of O_OX
    elif(ffatype == ffatype_3): ffidx_3 = idx     # determine index of O_OH
    if(idx > ffidx_max): ffidx_max = idx          # determine maximum index
    print(idx, ffatype)

print(ffatype_ids_copy)
print(ffidx_1, ffidx_2, ffidx_3)

# Determine the indices of the atoms belonging to one of the three atom types that will be changed
all_atoms = np.arange(len(system.numbers))
atoms_ffatype_1 = all_atoms[system.ffatype_ids == ffidx_1]
atoms_ffatype_2 = all_atoms[system.ffatype_ids == ffidx_2]
atoms_ffatype_3 = all_atoms[system.ffatype_ids == ffidx_3]

# Renumber the ffatype_ids, so that all indices except for the ones we are going to change are already correct
ffatype_ids_copy -= 1*(ffatype_ids_copy > ffidx_1) + 1*(ffatype_ids_copy > ffidx_2) + 1*(ffatype_ids_copy > ffidx_3)
ffidx_max -= 3

# Make the new ffatypes object
ffatypes_new = [ffatype for ffatype in system.ffatypes if ffatype not in (ffatype_1, ffatype_2, ffatype_3)]
ffatypes_new.extend(('O_CA_ZR', 'O_CA_HF', 'O_OX_0HF', 'O_OX_1HF', 'O_OX_2HF', 'O_OX_3HF', 'O_OH_0HF', 'O_OH_1HF', 'O_OH_2HF', 'O_OH_3HF'))

# Split the O_CA, depending on whether they are bonded to a zirconium or a cerium/hafnium
for idx in atoms_ffatype_1:
    neighs = system.neighs1[idx]
    for idx_2 in neighs:
        if system.numbers[idx_2] == idx_zr:
            ffatype_ids_copy[idx] = ffidx_max+1
        elif system.numbers[idx_2] == idx_ce:
            ffatype_ids_copy[idx] = ffidx_max+2

# Split the O_OX, depending on the number of neighbouring cerium/hafnium atoms
for idx in atoms_ffatype_2:
    neighs = system.neighs1[idx]
    counter = 0
    for idx_2 in neighs:
        if system.numbers[idx_2] == idx_ce:
            counter += 1
    if counter == 0: ffatype_ids_copy[idx] = ffidx_max+3
    elif counter == 1: ffatype_ids_copy[idx] = ffidx_max+4
    elif counter == 2: ffatype_ids_copy[idx] = ffidx_max+5
    elif counter == 3: ffatype_ids_copy[idx] = ffidx_max+6

# Split the O_OH, depending on the number of neighbouring cerium/hafnium atoms
for idx in atoms_ffatype_3:
    neighs = system.neighs1[idx]
    counter = 0
    for idx_2 in neighs:
        if system.numbers[idx_2] == idx_ce:
            counter += 1
    if counter == 0: ffatype_ids_copy[idx] = ffidx_max+7
    elif counter == 1: ffatype_ids_copy[idx] = ffidx_max+8
    elif counter == 2: ffatype_ids_copy[idx] = ffidx_max+9
    elif counter == 3: ffatype_ids_copy[idx] = ffidx_max+10

print(ffatype_ids_copy)

ff_idx_to_be_deleted = []

for ff_idx, ffatype in enumerate(ffatypes_new):
    # Check whether this atom type actually appears in the system
    if sum(ffatype_ids_copy == ff_idx) == 0:
        ff_idx_to_be_deleted.append(ff_idx)
        print('Will delete %s, as this atom type does not appear in the system' %(ffatype))

# Delete all unnecessary atom types
for ff_idx in sorted(ff_idx_to_be_deleted, reverse=True):
    # Delete from the ff atom type object
    ffatypes_new.pop(ff_idx)
    # Subtract one from subsequent ff atom type indices
    ffatype_ids_copy -= 1*(ffatype_ids_copy > ff_idx)

# Write out the new chk file
system.ffatypes = ffatypes_new
system.ffatype_ids = ffatype_ids_copy
system.to_file('init_new.chk')

