#!/usr/bin/env python
# File name: comparison_3x3x3.py
# Author: Joachim Vandewalle
# Date: 01-10-2022

"""Create a bar chart of the atomic and micromechanical bulk modulus for 3x3x3 systems."""

import json

import numpy as np
import pickle as pkl
import matplotlib.pyplot as plt


with open("output_configurations.json") as jfile:
    jdict = json.load(jfile)

x = np.arange(10)
confs = [str(i) for i in range(10)]
plt.rcParams['axes.axisbelow'] = True
#plt.grid()
plt.plot(x, [t[0] for t in jdict["all_E_atomic"]], label="E_min (atomic)")
plt.plot(x, [t[1] for t in jdict["all_E_atomic"]], label="E_max (atomic)")
plt.plot(x, [t[0] for t in jdict["all_E_micmec"]], label="E_min (micmec)")
plt.plot(x, [t[1] for t in jdict["all_E_micmec"]], label="E_max (micmec)")
plt.xticks(x, confs)
plt.xlim(-0.60, 9.60)
plt.ylim(0.0, 60.0)
plt.legend()

plt.axvspan(1.5, 4.5, facecolor="lightgrey", alpha=0.5, zorder=-100)
plt.axvspan(5.5, 8.5, facecolor="lightgrey", alpha=0.5, zorder=-100)
plt.ylabel("YOUNG MODULUS [GPa]")
plt.xlabel("CONFIGURATIONS")
plt.savefig("fig_young_modulus.png")
plt.show()

