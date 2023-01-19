#!/usr/bin/env python
# File name: fig_timesteps_2x2x2_fcu.py
# Author: Joachim Vandewalle
# Date: 20-11-2021

import numpy as np
import h5py

import matplotlib.pyplot as plt

from molmod.units import femtosecond, picosecond, kjmol

timesteps = np.linspace(100.0, 500.0, 5)

colors = [
    "#e41a1c", "#377eb8", "#4daf4a", 
    "#984ea3", "#ff7f00", "#ffff33", "#a65628",
    "#f781bf", "#999999", "#66c2a5", "#fc8d62",
    "#8da0cb", "#e78ac3", "#a6d854", "#ffd92f",
    "#e5c494", "#b3b3b3", "#8dd3c7", "#ffffb3",
    "#bebada", "#fb8072", "#80b1d3", "#fdb462",
    "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd",
    "#ccebc5", "#ffed6f"
]

for idx, timestep in enumerate(timesteps[:]):
    
    h5_fn = f"timestep{idx}.h5"

    with h5py.File(h5_fn, mode = 'r') as h5_f:
        
        dset_traj_energy = h5_f['trajectory/etot']
        dset_traj_energy_pot = h5_f['trajectory/epot']
        dset_traj_energy_kin = h5_f['trajectory/ekin']
        dset_traj_cons_err = h5_f['trajectory/cons_err']
                
        # Copy datasets to arrays.
        traj_energy = np.array(dset_traj_energy)
        traj_energy_pot = np.array(dset_traj_energy_pot)
        traj_energy_kin = np.array(dset_traj_energy_kin)
        traj_cons_err = np.array(dset_traj_cons_err)

        time = timestep*np.arange(len(traj_energy))*femtosecond
        
        plt.plot(time/picosecond, traj_energy/kjmol, label=f"{int(timestep)/1000} ps", color=colors[idx], linestyle="-")

plt.xlabel("TIME [ps]")
plt.ylabel("ENERGY [kJ/mol]")
plt.xlim(0.0, 10.0)
plt.ylim(0.0, 100.0)
plt.legend()
plt.grid()
plt.savefig("fig_timesteps_2x2x2_test.svg", format="svg")
plt.show()

