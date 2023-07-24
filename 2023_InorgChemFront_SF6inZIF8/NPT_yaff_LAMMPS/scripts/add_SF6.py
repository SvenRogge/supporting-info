from yaff import System
import numpy as np
from molmod.units import angstrom

sys_ZIF8_unit = System.from_file('ZIF8_empty.chk')
sys_ZIF8_empty = sys_ZIF8_unit.supercell(2,2,2)
n_atoms_ZIF8 = len(sys_ZIF8_empty.masses)
n_ffatypes_ZIF8 = len(sys_ZIF8_empty.ffatypes)

sys_SF6 = System.from_file('SF6.chk')
n_atoms_SF6 = len(sys_SF6.masses)
bonds_SF6 = sys_SF6.bonds.copy()

rvecs = sys_ZIF8_unit.cell.rvecs

origin = 0.5*np.diag(rvecs)

pos_SF6_centered = np.array(sys_SF6.pos) - np.array(sys_SF6.pos)[0,:]


displacements = np.array([[-3.5,0,0],[3.5,0,0],[0,3.5,0],[0,-3.5,0],[0,0,3.5],[0,0,-3.5]])*angstrom
n_mols = 6

pos_all = sys_ZIF8_empty.pos.copy()
bonds_all = sys_ZIF8_empty.bonds.copy()

for i in range(n_mols):
    disp = displacements[i]
    pos_all = np.vstack((pos_all, pos_SF6_centered + disp + origin))
    bonds_all = np.vstack((bonds_all, bonds_SF6 + n_atoms_ZIF8 + i* n_atoms_SF6))
    
    
def concatenate_same(array_1, array_2, n2, offset2=0):
    return np.concatenate((array_1, np.tile(array_2 + offset2, n2)))

ffatypes_all = list(sys_ZIF8_empty.ffatypes) + list(sys_SF6.ffatypes)
charges_all = concatenate_same(sys_ZIF8_empty.charges, sys_SF6.charges, n_mols)
masses_all = concatenate_same(sys_ZIF8_empty.masses, sys_SF6.masses, n_mols)
numbers_all = concatenate_same(sys_ZIF8_empty.numbers, sys_SF6.numbers, n_mols)
radii_all = concatenate_same(sys_ZIF8_empty.radii, sys_SF6.radii, n_mols)
ffatype_ids_all = concatenate_same(sys_ZIF8_empty.ffatype_ids, sys_SF6.ffatype_ids, n_mols, offset2 = n_ffatypes_ZIF8)

sys_filled = System(numbers=numbers_all, pos=pos_all, ffatypes=ffatypes_all, ffatype_ids=ffatype_ids_all, bonds=bonds_all, rvecs=sys_ZIF8_empty.cell.rvecs, charges=charges_all, radii=radii_all, masses=masses_all)
sys_filled.to_file('ZIF8_%dSF6.chk' %n_mols)

print(np.mean(sys_ZIF8_empty.pos, axis=0))
print(origin)
