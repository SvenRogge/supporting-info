#!/usr/bin/env python
# File name: comparison_3x3x3.py
# Author: Joachim Vandewalle
# Date: 01-10-2022

"""Create a 3D plot of the PES of a micromechanical system."""

import json

import numpy as np
import pickle as pkl
import matplotlib.pyplot as plt


with open("output_energy.json") as jfile:
    jdict = json.load(jfile)

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, projection='3d')

ax.plot_surface(np.array(jdict["X_ana"]), np.array(jdict["Y_ana"]), np.array(jdict["ENERGY_POT"]), cmap="viridis")
ax.set_xlabel("$x - x_0$ [Å]")
ax.set_ylabel("$y - y_0$ [Å]")
ax.set_zlabel("POTENTIAL ENERGY [kJ/mol]")

plt.savefig("output_energy.png", format="png")
plt.show()

