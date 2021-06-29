import numpy as np
import h5py
from yaff import System, Cell
import matplotlib.pyplot as pt
from molmod.units import picosecond, angstrom

sixrings = [[223, 220, 226, 217, 224, 218],
            [217, 224, 219, 223, 220, 227],
            [216, 227, 220, 222, 219, 224],
            [216, 224, 218, 222, 220, 226],
            [217, 227, 221, 223, 219, 225],
            [226, 216, 225, 218, 222, 221],
            [223, 218, 225, 217, 226, 221],
            [219, 225, 216, 227, 221, 222]]

mids_ring_init = np.array([[0.25,0.25,-0.25],[-0.25,0.25,-0.25],[-0.25,0.25,0.25],[0.25,0.25,0.25],[-0.25,-0.25,-0.25],[0.25,-0.25,0.25],[0.25,-0.25,-0.25],[-0.25,-0.25,0.25]])

n_framework_atoms = 276

h5_fns = ['traj_%d.h5' %i for i in range(1,21)]

with h5py.File(h5_fns[0], 'r') as f:
    numbers = np.array(f['system/numbers'])
    masses = np.array(f['system/masses'])
    pos = np.array(f['trajectory/pos'])
    cell = np.array(f['trajectory/cell'])

for h5_fn in h5_fns[1:]:
    print('Processing %s' %h5_fn)
    with h5py.File(h5_fn,'r') as f:
        pos = np.concatenate((pos, np.array(f['trajectory/pos'])))
        cell = np.concatenate((cell, np.array(f['trajectory/cell'])))

n_water = int((len(numbers)-n_framework_atoms)/3)
print("Found %i water molecules" %(n_water))

avg_pos = pos.mean(axis=1)
print(avg_pos[0,:]/angstrom)

pos_sixrings = np.zeros((pos.shape[0], len(sixrings), 3))

for t in range(pos.shape[0]):
    cell_t = Cell(cell[t,:,:])
    for idx_r, sixring in enumerate(sixrings):
        pos_0 = np.dot(mids_ring_init[idx_r], cell_t.rvecs) + avg_pos[t,:]
        pos_sixring = pos_0[:]
        for zn_id in sixring:
            delta = pos[t,zn_id,:] - pos_0
            cell_t.mic(delta)
            pos_sixring += delta/6
        pos_sixrings[t,idx_r,:] = pos_sixring


if n_water != 0:
    O_indices = range(n_framework_atoms, pos.shape[1], 3)
    pos_new = np.concatenate((pos_sixrings, pos[:,O_indices,:]), axis=1)
    numbers_new = np.concatenate((84*np.ones(len(sixrings)), numbers[O_indices]))
    masses_new = np.concatenate((np.zeros(len(sixrings)), masses[O_indices]))
else:
    pos_new = pos_sixrings
    numbers_new = 84*np.ones(len(sixrings))
    masses_new = np.zeros(len(sixrings))

with h5py.File('rdf.h5', 'w') as g:
    sgrp_g = g.create_group('system')
    sgrp_g.create_dataset('numbers', data=numbers_new)
    sgrp_g.create_dataset('masses', data=masses_new)
    sgrp_g.create_dataset('pos', data=pos_new[0,:,:])
    sgrp_g.create_dataset('rvecs', data=cell[0,:,:])

    tgrp_g = g.create_group('trajectory')
    tgrp_g.create_dataset('pos', data=pos_new)
    tgrp_g.create_dataset('cell', data=cell)
