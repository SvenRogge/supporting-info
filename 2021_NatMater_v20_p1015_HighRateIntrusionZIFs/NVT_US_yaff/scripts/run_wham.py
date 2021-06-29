import numpy as np
from molmod.units import *
import matplotlib.pyplot as pl
import h5py, os

class concat_h5(object):

    def __init__(self, fn_prefix='traj', directory='.'):

        self.keys = ['system/numbers', 'trajectory/time', 'trajectory/pos', 'trajectory/cell', 'trajectory/cv_values']
        self.h5py = {}

        fns = []
        for fn in os.listdir(directory):
            if fn.startswith(fn_prefix) and fn.endswith('.h5'):
                fns.append(fn)

        fns = sorted(fns, key=lambda a: int(a.split('.')[0].split('_')[1]))

        f = h5py.File('{}/{}'.format(directory, fns[0]), mode='r')
        self.numbers = np.array(f['system/numbers'])
        self.pos = np.array(f['trajectory/pos'])
        self.cell = np.array(f['trajectory/cell'])
        self.time = np.array(f['trajectory/time'])
        self.cv_values = np.array(f['trajectory/cv_values'])

        for fn in fns[1:]:
            g = h5py.File('{}/{}'.format(directory, fn), 'r')
            self.pos = np.concatenate((self.pos, g['trajectory/pos'][1:]), axis=0)
            self.cell = np.concatenate((self.cell, g['trajectory/cell'][1:]), axis=0)
            self.time = np.concatenate((self.time, g['trajectory/time'][1:]), axis=0)
            self.cv_values = np.concatenate((self.cv_values, g['trajectory/cv_values'][1:]), axis=0)

        for key in self.keys:
            self.h5py[key] = eval('self.{}'.format(key.split('/')[1]))

    def __setitem__(self, key, data):
        self.h5py[key] = data

    def __getitem__(self, key):
        return self.h5py[key]


# Generate WHAM input files
cv_start = -11.5
cv_stop = 11.5
cv_step = 0.5
eq_step = 200

kappa = 25 # kjmol / [CV]^2 (WHAM set to kJ/mol)

f = open('metadata.dat', 'w')
f.write('# path_timeseries_file\tCV0\tkappa\n')

for cv in np.arange(cv_start, cv_stop+0.5*cv_step, cv_step):

    dirname = '{:.1f}'.format(cv)
    if cv < 0: dirname = 'm' + dirname[1:]

    if not os.path.exists('colvar_{}'.format(dirname)):

        g = concat_h5(directory=dirname)
        time = g['trajectory/time'][eq_step:]/picosecond
        cv_values = g['trajectory/cv_values'][eq_step:]/angstrom

        np.savetxt('colvar_{}'.format(dirname), np.c_[time, cv_values], header='Time (ps)\tCV (angstrom)')

    f.write('colvar_{}\t{}\t{}\n'.format(dirname, cv, kappa))

f.close()

# Execute WHAM
temp = 300
bin_width = 0.05
error_analysis = False

if not error_analysis:
    os.system('wham {} {} {} 1e-6 {} 0 metadata.dat Free_energy.dat'.format(cv_start, cv_stop, int((cv_stop-cv_start)/bin_width), temp))

else:
    MC_trials = 50
    os.system('wham {} {} {} 1e-6 {} 0 metadata.dat Free_energy.dat {} 12321'.format(cv_start, cv_stop, int((cv_stop-cv_start)/bin_width), temp, MC_trials))

# Plot free energy
F = np.loadtxt('Free_energy.dat')
pl.plot(F[:, 0], F[:, 1])
pl.ylim([-5, 30])
pl.xlabel('CV (angstrom)')
pl.ylabel('F (kJ/mol)')
pl.savefig('F.pdf', bbox_inches='tight')
pl.show()
