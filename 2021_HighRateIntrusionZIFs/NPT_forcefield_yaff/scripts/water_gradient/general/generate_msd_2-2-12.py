import numpy as np
import h5py
from yaff import System, Cell
import matplotlib.pyplot as pt
from molmod.units import picosecond, angstrom, nanosecond


n_z = 12
n_traj = 10

frac_coords_pores_basis = np.array([[0.0,0.0,0.0],[0.5,0.0,0.0],[0.0,0.5,0.0],[0.5,0.5,0.0],[0.25,0.25,1./(2*n_z)],[0.75,0.25,1./(2*n_z)],[0.25,0.75,1./(2*n_z)],[0.75,0.75,1./(2*n_z)]])
offsets = np.array([[0,0,1.0*x/n_z] for x in range(12)])
frac_coords_pores = np.array([c+offset for offset in offsets for c in frac_coords_pores_basis])

n_framework_atoms = 2*2*n_z*276

with h5py.File('traj_1.h5', 'r') as f:
    numbers = np.array(f['system/numbers'])

n_water = int((len(numbers)-n_framework_atoms)/3)

## Find the water molecules
O_water = []
water_idx = []

for idx_n, number in enumerate(numbers):
    if number == 8 and idx_n < len(numbers)-2:
	if numbers[idx_n+1] == 1 and numbers[idx_n+2] ==1:
		O_water.append(idx_n)
                water_idx.append(idx_n)
		water_idx.append(idx_n+1)
		water_idx.append(idx_n+2)

assert len(O_water) == n_water
print("Found %i water molecules" %(len(O_water)))

framework_idx = [x for x in range(len(numbers)) if x not in water_idx]
# Check that the water molecules are indeed stored as [O, H, H]
#for i in range(n_framework_atoms, len(numbers), 3):
#    assert numbers[i] == 8


def load_h5(h5_name, time, pos_molecules, cell, ref_origin, n_framework_atoms, numbers):
	with h5py.File(h5_name, 'r') as f:
		time = np.append(time, np.array(f['trajectory/time']))
		pos_molecules = np.append(pos_molecules, np.array(f['trajectory/pos'][:,O_water,:]), axis=0)
		cell = np.append(cell, np.array(f['trajectory/cell']), axis=0)
                ref_origin = np.append(ref_origin, np.mean(np.array(f['trajectory/pos'][:,framework_idx,:]),axis=1), axis=0)
	return time, pos_molecules, cell, ref_origin

with h5py.File('traj_1.h5', 'r') as f:
	time = np.array(f['trajectory/time'])
	pos_molecules = np.array(f['trajectory/pos'][:,O_water,:])
	cell = np.array(f['trajectory/cell'])
        ref_origin = np.mean(np.array(f['trajectory/pos'][:,framework_idx,:]),axis=1)
        print(ref_origin.shape)

for i in range(1,n_traj):
	time, pos_molecules, cell, ref_origin = load_h5('traj_%s.h5' %(i+1), time, pos_molecules, cell, ref_origin, n_framework_atoms, numbers)

print(ref_origin.shape)

print(pos_molecules.shape)
loc_molecule = np.zeros((pos_molecules.shape[0], pos_molecules.shape[1]))
pore_pos = np.zeros((pos_molecules.shape[0], frac_coords_pores.shape[0], 3))
min_dist_all = np.zeros((pos_molecules.shape[0], pos_molecules.shape[1]))

for t in range(pos_molecules.shape[0]):
    cell_t = Cell(cell[t,:,:])
    cutoff = np.sqrt(3)*7*angstrom
    for idx_p, pore_frac in enumerate(frac_coords_pores):
        pore_pos[t,idx_p,:] = np.dot(pore_frac, cell_t.rvecs)+ref_origin[t,:]-ref_origin[0,:]

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
        min_dist_all[t,idx_m] = min_dist

with h5py.File('msd_long.h5', 'w') as f:
    f.create_dataset('time', data=time)
    f.create_dataset('pos_molecules', data=pos_molecules)
    f.create_dataset('pos_pores', data=pore_pos)
    f.create_dataset('loc_molecules', data=loc_molecule)
    f.create_dataset('min_dist', data=min_dist_all)
    f.create_dataset('ref_origin', data=ref_origin)


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

