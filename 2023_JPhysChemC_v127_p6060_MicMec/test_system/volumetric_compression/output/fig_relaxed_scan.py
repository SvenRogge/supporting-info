#!/usr/bin/env python
# File name: comparison_3x3x3.py
# Author: Joachim Vandewalle
# Date: 01-10-2022

"""Plot the volume-energy relaxed scan of a 2x2x2 and 3x3x3 test system."""

import json

import numpy as np
import pickle as pkl
import matplotlib.pyplot as plt


with open("output_relaxed_scan.json") as jfile:
    jdict = json.load(jfile)

idx = 0
for scalings, energy_density in zip(jdict["all_scalings"], jdict["all_energy_densities"]):
    if idx == 0:
        plt.plot(scalings, energy_density, label="2 x 2 x 2")
    else:
        plt.plot(scalings, energy_density, "--", label="3 x 3 x 3")
    idx += 1
    
#plt.plot(scalings, 0.5**(13*gigapascal/np.linalg.det(orig_rvecs))*(scalings - 1.0)**2)
plt.xlabel(r"$V/V_0$ [-]")
plt.ylabel("POTENTIAL ENERGY DENSITY [kJ/mol/Å³]")
plt.xlim(scalings[0], scalings[-1])
plt.ylim(0.0)
plt.legend()
plt.grid()
plt.savefig("fig_relaxed_scan.svg", format="svg")
plt.show()

