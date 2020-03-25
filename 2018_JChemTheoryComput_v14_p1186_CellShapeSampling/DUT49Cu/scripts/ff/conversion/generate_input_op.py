#! /usr/bin/env python

from molmod.units import *
from yaff import *
import numpy as np
import h5py
from molmod.periodic import periodic

f_input = "DUT-49op-DFT.xyz"
f_output = "output.arc"

def determine_type(atom_number, neighs, neighs2):
    if atom_number == 1:
        return 5    # hydrogen
    elif atom_number == 7:
        return 40   # nitrogen
    elif atom_number == 8:
        return 167  # oxygen
    elif atom_number == 29:
        return 165  # copper
    else: # carbon
        neighs_num = [sys.numbers[i] for i in neighs]
        if 8 in neighs_num: 
            return 168  # carboxylate carbon
        elif 1 in neighs_num or 7 in neighs_num:
            return 2    # regular carbon
        else:
            neighs2_num = [sys.numbers[i] for i in neighs2]
            if 7 in neighs2_num or 8 in neighs2_num:
                return 2    # regular carbon
            else:
                return 180  # bent carbon

# Create system object
sys = System.from_file(f_input,rvecs=46.72*angstrom*np.eye(3))
sys.detect_bonds()
print sys.cell.rvecs
# sys.neighs now contains first neighbors, sys.neighs2 contains second neighbors
# sys.numbers contains the atom numbers
# sys.pos contains the positions

# Initialize list with atom type
atypes = []

# Initialize list with neighbors
nlist = []

for i in xrange(len(sys.numbers)):
    new_type = determine_type(sys.numbers[i], sys.neighs1[i], sys.neighs2[i])
    atypes.append(new_type)
    nlist.append(list(sys.neighs1[i]))

# We still need to add "neighboring" coppers. So run over carboxylate carbons
for i in xrange(len(sys.numbers)):
    if atypes[i] == 168:
        neighs2 = sys.neighs2[i]
        pair = []
        for idx in neighs2:
            if sys.numbers[idx] == 29:
                pair.append(idx)
        # pair should contain 2 coppers now
        if not pair[1] in nlist[pair[0]]:
            nlist[pair[0]].append(pair[1])
        if not pair[0] in nlist[pair[1]]:
            nlist[pair[1]].append(pair[0])

for i in xrange(len(sys.numbers)):
    if periodic[sys.numbers[i]].symbol == "H":
        if len(nlist[i]) != 1:
            print("%d %s" % (i+1, periodic[sys.numbers[i]].symbol))
    elif periodic[sys.numbers[i]].symbol == "C":
        if len(nlist[i]) != 3:
            print("%d %s" % (i+1, periodic[sys.numbers[i]].symbol))
    elif periodic[sys.numbers[i]].symbol == "O":
        if len(nlist[i]) != 2:
            print("%d %s" % (i+1, periodic[sys.numbers[i]].symbol))
    elif periodic[sys.numbers[i]].symbol == "Cu":
        if len(nlist[i]) != 5:
            print("%d %s" % (i+1, periodic[sys.numbers[i]].symbol))
    

# Write out to xyz
f = open(f_output, 'w')
f.write("1728   46.723700000000001  46.723700000000001  46.723700000000001    90.0  90.0  90.0 \n")
for i in xrange(len(sys.numbers)):
    pos = sys.pos[i]
    string_begin = "%d %s\t%.8f\t%.8f\t%.8f\t%d" %(i+1, periodic[sys.numbers[i]].symbol, pos[0]/angstrom, pos[1]/angstrom, pos[2]/angstrom, atypes[i])
    for item in nlist[i]:
        string_begin = string_begin + "\t%d" %item
    f.write(string_begin + "\n")

f.close()
