import numpy as np
import h5py
from yaff import System, Cell
import matplotlib.pyplot as pt
from molmod.units import picosecond, angstrom, nanosecond


frac_coords_pores = np.array([[0.0,0.0,0.0],[0.5,0.5,0.5]])
n_framework_atoms = 276


with h5py.File('msd_long.h5', 'r') as f:
    loc_molecule = np.array(f['loc_molecules'])
    pos_molecule = np.array(f['pos_molecules'])
    time = np.array(f['time'])

for i in range(pos_molecule.shape[1]):
    pos_molecule[:,i,:] -= pos_molecule[0,i,:]

n_molecules = pos_molecule.shape[1]

msd_all = np.zeros((pos_molecule.shape[0], pos_molecule.shape[1]))

print(pos_molecule.shape)

for t in range(pos_molecule.shape[0]):
    dt = t+1
    # Run over every molecule
    for j in range(pos_molecule.shape[1]):
        msd_tot = 0
        counter = 0
        # Run over every time interval
        for i in range(len(time)//dt):
            msd_tot += ((pos_molecule[(i+1)*dt-1,j,:] - pos_molecule[i*dt,j,:])**2).sum()
            counter += 1
        msd_all[t,j] = msd_tot/counter

pt.clf()
pt.xlabel('Simulation time [ns]')
pt.ylabel('Mean squared displacement [\AA$^2$]')

pt.plot(time/nanosecond,(pos_molecule**2).sum(axis=1).sum(axis=1)/n_molecules/(angstrom**2), '0.5',label='total')
pt.plot(time/nanosecond,(pos_molecule[:,:,0]**2).sum(axis=1)/n_molecules/(angstrom**2), '#ef6548', label='x component')
pt.plot(time/nanosecond,(pos_molecule[:,:,1]**2).sum(axis=1)/n_molecules/(angstrom**2), '#41ab5d', label='y component')
pt.plot(time/nanosecond,(pos_molecule[:,:,2]**2).sum(axis=1)/n_molecules/(angstrom**2), '#3690c0', label='z component')
pt.plot(time[:len(time)/2]/nanosecond,msd_all[:len(time)/2,:].sum(axis=1)/n_molecules/(angstrom**2), 'k', label='multiple origin')

legend = pt.legend(loc='upper left',fancybox=True)
pt.xlim([0,5])
pt.tight_layout()
pt.savefig('ZIF8_diff_H2O_long.svg', format='svg', bbox_inches = 'tight')
pt.savefig('ZIF8_diff_H2O_long.png', bbox_inches = 'tight')


pt.clf()
pt.xlabel('Simulation time [ns]')
pt.ylabel('Diffusion coefficient [\AA$^2$/nanosecond]')

pt.plot(time[1000:]/nanosecond,(pos_molecule[1000:,:,:]**2).sum(axis=1).sum(axis=1)/(6*time[1000:])/n_molecules/(angstrom**2/nanosecond), '0.5',label='total')
pt.plot(time[1000:]/nanosecond,(pos_molecule[1000:,:,0]**2).sum(axis=1)/(6*time[1000:])/n_molecules/(angstrom**2/nanosecond), '#ef6548', label='x component')
pt.plot(time[1000:]/nanosecond,(pos_molecule[1000:,:,1]**2).sum(axis=1)/(6*time[1000:])/n_molecules/(angstrom**2/nanosecond), '#41ab5d', label='y component')
pt.plot(time[1000:]/nanosecond,(pos_molecule[1000:,:,2]**2).sum(axis=1)/(6*time[1000:])/n_molecules/(angstrom**2/nanosecond), '#3690c0', label='z component')
pt.plot(time[1000:len(time)/2]/nanosecond,(msd_all[1000:len(time)/2,:]).sum(axis=1)/n_molecules/(6*time[1000:len(time)/2])/(angstrom**2/nanosecond), 'k', label='multiple origin')

legend = pt.legend(loc='upper right',fancybox=True)
pt.xlim([0,5])
pt.tight_layout()
pt.savefig('ZIF8_D_H2O_long.svg', format='svg', bbox_inches = 'tight')
pt.savefig('ZIF8_D_H2O_long.png', bbox_inches = 'tight')

