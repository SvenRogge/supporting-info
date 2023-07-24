import numpy as np
from molmod.units import angstrom, picosecond
import matplotlib.pyplot as pt
import h5py, os


# Generate WHAM input files
cv_start = -11
cv_stop = 11
cv_step = 0.5
eq_step = 200

kappas = [25, 50, 75, 100, 150, 200, 300] # kjmol / [CV]^2 (WHAM set to kJ/mol)

with open('../results/metadata.dat', 'w') as f:
    f.write('# path_timeseries_file\tCV0\tkappa\n')

    for kappa in kappas:
        if kappa > 25:
            cv_start_it = -3 # reduced cv sweep
            cv_stop_it = 3 # reduced cv sweep
        else:
            cv_start_it = cv_start
            cv_stop_it = cv_stop
        for cv in np.arange(cv_start_it, cv_stop_it+0.5*cv_step, cv_step):

            dirname = '%d' %(10*cv)
            if cv < 0: dirname = 'm' + dirname[1:]
            if kappa > 25: dirname = dirname + '_%d' %kappa
            print(dirname)

            if not os.path.exists('../results/colvar_%s' %dirname):
                with h5py.File('../%s/traj_1.h5' %dirname, 'r') as g:
                    time = g['trajectory/time'][eq_step:]
                    cv_values = g['trajectory/cv_values'][eq_step:]
                np.savetxt('../results/colvar_%s' %dirname, np.c_[time/picosecond, cv_values/angstrom], header='Time (ps)\tCV (angstrom)')

            f.write('../results/colvar_%s\t%.1f\t%d\n' %(dirname, cv, kappa))


# Execute WHAM
temp = 300
bin_width = 0.01
error_analysis = False

if not error_analysis:
    os.system('wham {} {} {} 1e-6 {} 0 ../results/metadata.dat ../results/free_energy.dat'.format(cv_start, cv_stop, int((cv_stop-cv_start)/bin_width), temp))

else:
    MC_trials = 50
    os.system('wham {} {} {} 1e-6 {} 0 ../results/metadata.dat ../results/free_energy.dat {} 12321'.format(cv_start, cv_stop, int((cv_stop-cv_start)/bin_width), temp, MC_trials))

# Plot free energy
free_energy = np.loadtxt('../results/free_energy.dat')
pt.plot(free_energy[:, 0], free_energy[:, 1])
pt.xlabel('CV (Angstrom)')
pt.ylabel('F (kJ/mol)')
pt.savefig('../results/free_energy.svg', bbox_inches='tight')
pt.show()
