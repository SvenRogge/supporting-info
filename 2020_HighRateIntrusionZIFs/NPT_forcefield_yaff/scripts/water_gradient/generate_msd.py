import numpy as np
import h5py
from yaff import System, Cell
import matplotlib.pyplot as pt
from molmod.units import picosecond, angstrom, nanosecond

frac_coords_pores = np.array([[0.0,0.0,0.0],[0.5,0.5,0.25],[0.0,0.0,0.5],[0.5,0.5,0.75]])
n_framework_atoms = 2*276

with h5py.File('traj_1.h5', 'r') as f:
    numbers = np.array(f['system/numbers'])

n_water = int((len(numbers)-n_framework_atoms)/3)

print("Found %i water molecules" %(n_water))

# Check that the water molecules are indeed stored as [O, H, H]
#for i in range(n_framework_atoms, len(numbers), 3):
#    assert numbers[i] == 8


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
    cutoff = 0.4*np.mean([cell[t,i,i] for i in range(2)])
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


'''

with h5py.File('msd.h5', 'r') as f:
    loc_molecule = np.array(f['loc_molecules'])
    pos_molecule = np.array(f['pos_molecules'])
    time = np.array(f['time'])

# Identify possible hopping
hopping = False
for idx, molecule_loc in enumerate(loc_molecule.T):
    old_value = molecule_loc[0]
    for idx2, new_value in enumerate(molecule_loc):
        if old_value != new_value and new_value != -1:
            old_value = new_value
            hopping = True
            print('Hopping of molecule %i (O index: %i) at time step %i' %(idx, n_framework_atoms+3*idx, idx2))
if not hopping: print('No hopping event encountered')

for i in range(pos_molecules.shape[1]):
    pos_molecules[:,i,:] -= pos_molecules[0,i,:]

msd_all = np.zeros((pos_molecules.shape[0], pos_molecules.shape[1]))

print(pos_molecules.shape)

for t in range(pos_molecules.shape[0]):
    dt = t+1
    # Run over every molecule
    for j in range(pos_molecule.shape[1]):
        msd_tot = 0
        counter = 0
        # Run over every time interval
        for i in range(len(time)//dt):
            msd_tot += ((pos_molecules[(i+1)*dt-1,j,:] - pos_molecules[i*dt,j,:])**2).sum()
            counter += 1
        msd_all[i,j] = msd_tot/counter

pt.clf()
pt.rc('text', usetex=True)
pt.rc('font',**{'family':'sans-serif','sans-serif':['Paladino']})
pt.rcParams['font.family'] = 'Paladino'
comap = pt.cm.get_cmap(name='jet')
pt.xlabel(r'Simulation time [ps]')
pt.ylabel(r'MSD [\AA$^2$]')

pt.plot(time/picosecond,(pos_molecules**2).sum(axis=1).sum(axis=1)/(angstrom**2), '0.5',label='total')
pt.plot(time/picosecond,(pos_molecules[:,:,0]**2).sum(axis=1)/(angstrom**2), 'r', label='x component')
pt.plot(time/picosecond,(pos_molecules[:,:,1]**2).sum(axis=1)/(angstrom**2), 'g', label='y component')
pt.plot(time/picosecond,(pos_molecules[:,:,2]**2).sum(axis=1)/(angstrom**2), 'b', label='z component')
pt.plot(time/nanosecond,msd_all.sum(axis=1)/(angstrom**2), 'k', label='multiple origin')
#pt.xlim(0,15)

legend = pt.legend(loc='upper left',fancybox=True)
ltext = legend.get_texts()
pt.setp(ltext, fontsize='small')
pt.savefig('ZIF8_diff_H2O.pdf', format='pdf', bbox_inches = 'tight')
pt.savefig('ZIF8_diff_H2O.png', bbox_inches = 'tight')


pt.clf()
pt.rc('text', usetex=True)
pt.rc('font',**{'family':'sans-serif','sans-serif':['Paladino']})
pt.rcParams['font.family'] = 'Paladino'
comap = pt.cm.get_cmap(name='jet')
pt.xlabel(r'Simulation time [ps]')
pt.ylabel(r'D [\AA$^2$/nanosecond]')

pt.plot(time[10:]/picosecond,(pos_molecules[10:,:,:]**2).sum(axis=1).sum(axis=1)/(6*time[10:])/(angstrom**2/nanosecond), '0.5',label='total')
pt.plot(time[10:]/picosecond,(pos_molecules[10:,:,0]**2).sum(axis=1)/(6*time[10:])/(angstrom**2/nanosecond), 'r', label='x component')
pt.plot(time[10:]/picosecond,(pos_molecules[10:,:,1]**2).sum(axis=1)/(6*time[10:])/(angstrom**2/nanosecond), 'g', label='y component')
pt.plot(time[10:]/picosecond,(pos_molecules[10:,:,2]**2).sum(axis=1)/(6*time[10:])/(angstrom**2/nanosecond), 'b', label='z component')
pt.plot(time[10:]/picosecond,(msd_all[10:,:]).sum(axis=1)/(6*time[10:])/(angstrom**2/nanosecond), 'k', label='multiple origin')
#pt.xlim(0,15)

legend = pt.legend(loc='upper left',fancybox=True)
ltext = legend.get_texts()
pt.setp(ltext, fontsize='small')
pt.savefig('ZIF8_D_H2O.pdf', format='pdf', bbox_inches = 'tight')
pt.savefig('ZIF8_D_H2O.png', bbox_inches = 'tight')

'''

