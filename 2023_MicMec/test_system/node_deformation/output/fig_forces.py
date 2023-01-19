#!/usr/bin/env python
# File name: comparison_3x3x3.py
# Author: Joachim Vandewalle
# Date: 01-10-2022

"""Create a 3D plot of the PES of a micromechanical 2x2x2 fcu system."""

import json

import numpy as np
import pickle as pkl
import matplotlib.pyplot as plt


with open("output_forces.json") as jfile:
    jdict = json.load(jfile)

fig = plt.figure()

ax = fig.add_subplot(1, 1, 1)
ax.plot(np.array(jdict["X_ana"]), np.array(jdict["FX_ana"]), "o:", mfc="none", ms=16, label="analytical expression", color="blue")
ax.plot(np.array(jdict["X_num"]), np.array(jdict["FX_num"]), "x", ms=16, label="numerical derivative", color="red")
ax.set_xlabel("$x - x_0$ [Å]")
ax.set_ylabel("$f_x$ [kJ/mol/Å]")
ax.legend()
ax.grid()

plt.savefig("output_forces.png", format="png")
plt.show()

