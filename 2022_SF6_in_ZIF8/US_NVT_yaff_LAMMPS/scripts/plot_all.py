import numpy as np
import matplotlib.pyplot as pt


dirs = ['1SF6', '2SF6', '3SF6', '4SF6', '5SF6', '6SF6']

# Plot free energy
pt.clf()
for t_dir in dirs:  
    free_energy = np.loadtxt('%s/results/free_energy.dat' %t_dir)
    h = int(free_energy.shape[0]/2)
    pt.plot(free_energy[:, 0], free_energy[:, 1]-np.amin(free_energy[h:,1]), label=t_dir)
pt.xlabel('CV (Angstrom)')
pt.ylabel('F (kJ/mol)')
pt.legend()
pt.ylim([-40,120])
pt.savefig('free_energy_all.svg', bbox_inches='tight')
pt.show()
