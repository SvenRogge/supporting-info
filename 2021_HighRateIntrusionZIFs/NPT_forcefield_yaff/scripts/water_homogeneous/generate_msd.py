import numpy as np
import h5py
from yaff import System, Cell
import matplotlib.pyplot as pt
from molmod.units import picosecond, angstrom, nanosecond

frac_coords_pores = np.array([[0.0,0.0,0.0],[0.5,0.5,0.5]])
n_framework_atoms = 276

with h5py.File('traj_1.h5', 'r') as f:
    numbers = np.array(f['system/numbers'])

n_water = int((len(numbers)-n_framework_atoms)/3)

print("Found %i water molecules" %(n_water))

# Check that the water molecules are indeed stored as [O, H, H]
for i in range(n_framework_atoms, len(numbers), 3):
    assert numbers[i] == 8


def load_h5(h5_name, time, pos_molecules, cell, n_framework_atoms, numbers):
	with h5py.File(h5_name, 'r') as f:
		time = np.append(time, np.array(f['trajectory/time']))
		pos_molecules = np.append(pos_molecules, np.array(f['trajectory/pos'][:,n_framework_atoms:len(numbers):3,:]), axis=0)
		cell = np.append(cell, np.array(f['trajectory/cell']), axis=0)
	return time, pos_molecules, cell

with h5py.File('traj_1.h5', 'r') as f:
	time = np.array(f['trajectory/time'])
	pos_molecules = np.array(f['trajectory/pos'][:,n_framework_atoms:len(numbers):3,:])
	cell = np.array(f['trajectory/cell'])

for i in range(1,20):
	time, pos_molecules, cell = load_h5('traj_%s.h5' %(i+1), time, pos_molecules, cell, n_framework_atoms, numbers)

print(pos_molecules.shape)
loc_molecule = np.zeros((pos_molecules.shape[0], pos_molecules.shape[1]))
pore_pos = np.zeros((pos_molecules.shape[0], frac_coords_pores.shape[0], 3))

for t in range(pos_molecules.shape[0]):
    cell_t = Cell(cell[t,:,:])
    cutoff = 0.4*np.mean([cell[t,i,i] for i in range(3)])
    for idx_p, pore_frac in enumerate(frac_coords_pores):
        pore_pos[t,idx_p,:] = np.dot(pore_frac, cell_t.rvecs)

    for idx_m, pos_molecule in enumerate(pos_molecules[t,:,:]):
        min_dist = 100*angstrom
        loc_pore = -1
        for idx_p, pos_pore in enumerate(pore_pos[t,:,:]):
            rel_vec = pos_pore - pos_molecule
            cell_t.mic(rel_vec)
            dist = np.linalg.norm(rel_vec)
            if dist < min_dist and dist < cutoff:
                min_dist = dist
                loc_pore = idx_p
        loc_molecule[t,idx_m] = loc_pore

with h5py.File('msd_long.h5', 'w') as f:
    f.create_dataset('time', data=time)
    f.create_dataset('pos_molecules', data=pos_molecules)
    f.create_dataset('pos_pores', data=pore_pos)
    f.create_dataset('loc_molecules', data=loc_molecule)
