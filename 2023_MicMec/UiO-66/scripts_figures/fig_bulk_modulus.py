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
plt.bar(x - 0.20, jdict["all_K_atomic"], width=0.30)
plt.bar(x + 0.20, jdict["all_K_micmec"], width=0.30)
plt.xticks(x, confs)
plt.xlim(-0.60, 9.60)
plt.ylim(0.0, 40.0)
plt.legend(labels=["atomic", "micromechanical"])
# Bulk modulus of an fcu cell.
plt.hlines(jdict["K_fcu"], xmin=-0.60, xmax=9.60, linestyles="dashed", color="black")
# Bulk modulus of a reo cell.
plt.hlines(jdict["K_reo"], xmin=-0.60, xmax=9.60, linestyles="dotted", color="black")
plt.axvspan(1.5, 4.5, facecolor="lightgrey", alpha=0.5, zorder=-100)
plt.axvspan(5.5, 8.5, facecolor="lightgrey", alpha=0.5, zorder=-100)
plt.ylabel("BULK MODULUS [GPa]")
plt.xlabel("CONFIGURATIONS")
plt.savefig("fig_bulk_modulus.png")
plt.show()

