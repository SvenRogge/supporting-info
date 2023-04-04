#!/usr/bin/env python
# File name: compare_3x3x3_confs.py
# Author: Joachim Vandewalle
# Date: 01-10-2022

"""Compare elastic properties of 3x3x3 systems from the .pickle output of simulations."""

import json

import numpy as np
import pickle as pkl
import matplotlib.pyplot as plt

from micmec.analysis.tensor import *
from molmod.units import pascal
gigapascal = 1e9*pascal


fdir = ""

with open(fdir + "type_1x1x1_fcu_atomic.pickle", "rb") as f_fcu, open(fdir + "type_1x1x1_reo_atomic.pickle", "rb") as f_reo:
    data_fcu = pkl.load(f_fcu)
    data_reo = pkl.load(f_reo)
C_fcu = voigt(data_fcu["elasticity"][0], mode="elasticity")
C_reo = voigt(data_reo["elasticity"][0], mode="elasticity")
K_fcu = C_fcu[:3, :3].mean()
K_reo = C_reo[:3, :3].mean()

print("UiO-66(Zr): fcu")
print("---------------")
print("\nAtomic elasticity matrix and derived bulk modulus [GPa]:")
print(pretty_6x6_matrix(C_fcu/gigapascal), K_fcu/gigapascal)
print("\n")

print("UiO-66(Zr): reo")
print("---------------")
print("\nAtomic elasticity matrix and derived bulk modulus [GPa]:")
print(pretty_6x6_matrix(C_reo/gigapascal), K_reo/gigapascal)
print("\n")

confs = [str(i) for i in range(10)]
all_K_atomic = []
all_E_atomic = []
all_C_atomic = []

all_K_micmec = []
all_E_micmec = []
all_C_micmec = []


for conf in confs:

    print(f"CONFIGURATION{conf}")
    print("--------------")
    fn_atomic = f"type_3x3x3_conf{conf}_atomic.pickle"
    fn_micmec = f"type_3x3x3_conf{conf}_micmec.pickle"
    with open(fdir + fn_atomic, "rb") as f_atomic, open(fdir + fn_micmec, "rb") as f_micmec:
        data_atomic = pkl.load(f_atomic)
        data_micmec = pkl.load(f_micmec)
    C_atomic = voigt(data_atomic["elasticity"][0], mode="elasticity")
    C_micmec = voigt(data_micmec["elasticity"][0], mode="elasticity")
    
    all_C_micmec.append(C_micmec)
    all_C_atomic.append(C_atomic)

    K_atomic = bulk_modulus(C_atomic)
    all_K_atomic.append(K_atomic)
    all_E_atomic.append(min_max_young_modulus(C_atomic))
    
    K_micmec = bulk_modulus(C_micmec)
    all_K_micmec.append(K_micmec)
    all_E_micmec.append(min_max_young_modulus(C_micmec))

    print("\nAtomic elasticity matrix and derived bulk modulus [GPa]:")
    print(pretty_6x6_matrix(C_atomic/gigapascal), K_atomic/gigapascal)
    print("\nMicromechanical elasticity matrix and derived bulk modulus [GPa]:")
    print(pretty_6x6_matrix(C_micmec/gigapascal), K_micmec/gigapascal)

    # Elastic deformation modes.
    lC_atomic, vC_atomic = np.linalg.eig(C_atomic)
    idxs = lC_atomic.argsort()[::-1]
    lC_atomic = lC_atomic[idxs]
    vC_atomic = vC_atomic[:, idxs]

    lC_micmec, vC_micmec = np.linalg.eig(C_micmec)
    idxs = lC_micmec.argsort()[::-1]
    lC_micmec = lC_micmec[idxs]
    vC_micmec = vC_micmec[:, idxs]

    # Principal stresses and planes.
    for i in range(6):
        print(f"\nEigenvalue #{i} of atomic elasticity matrix [GPa]:")
        print(lC_atomic[i]/gigapascal)
        print(f"\nEigenvalue #{i} of micromechanical elasticity matrix [GPa]:")
        print(lC_micmec[i]/gigapascal)
        eps_atomic = voigt_inv(vC_atomic[:, i], mode="strain")
        eps_micmec = voigt_inv(vC_micmec[:, i], mode="strain")
        print(f"\nEigenvector #{i} of atomic elasticity matrix (strain tensor / elastic deformation mode):")
        print(eps_atomic)
        print(f"\nEigenvector #{i} of micromechanical elasticity matrix (strain tensor / elastic deformation mode):")
        print(eps_micmec)
        sigma_atomic = lC_atomic[i]*eps_atomic
        sigma_micmec = lC_micmec[i]*eps_micmec

        lsigma_atomic, vsigma_atomic = np.linalg.eig(sigma_atomic)
        idxs = lsigma_atomic.argsort()[::-1]
        lsigma_atomic = lsigma_atomic[idxs]
        vsigma_atomic = vsigma_atomic[:, idxs]

        lsigma_micmec, vsigma_micmec = np.linalg.eig(sigma_micmec)
        idxs = lsigma_micmec.argsort()[::-1]
        lsigma_micmec = lsigma_micmec[idxs]
        vsigma_micmec = vsigma_micmec[:, idxs]


        print("\nPrincipal stresses (atomic) [GPa]:")
        print(f"{lsigma_atomic[0]/gigapascal} ; {lsigma_atomic[1]/gigapascal} ; {lsigma_atomic[2]/gigapascal}")
        print("\nPrincipal planes (atomic):")
        print("[[{:6.3f}],   [[{:6.3f}],    [[{:6.3f}],".format(*list(vsigma_atomic[0])))
        print(" [{:6.3f}], ;  [{:6.3f}], ;   [{:6.3f}],".format(*list(vsigma_atomic[1])))
        print(" [{:6.3f}]]    [{:6.3f}]]     [{:6.3f}]]".format(*list(vsigma_atomic[2])))
        print("\nPrincipal stresses (micromechanical) [GPa]:")
        print(f"{lsigma_micmec[0]/gigapascal} ; {lsigma_micmec[1]/gigapascal} ; {lsigma_micmec[2]/gigapascal}")
        print("\nPrincipal planes (micromechanical):")
        print("[[{:6.3f}],   [[{:6.3f}],    [[{:6.3f}],".format(*list(vsigma_micmec[0])))
        print(" [{:6.3f}], ;  [{:6.3f}], ;   [{:6.3f}],".format(*list(vsigma_micmec[1])))
        print(" [{:6.3f}]]    [{:6.3f}]]     [{:6.3f}]]".format(*list(vsigma_micmec[2])))
        print("\n")
    print("\n")
    
    
C_micmec_conf0 = (all_C_micmec[0]/gigapascal).tolist()
C_micmec_conf2 = (all_C_micmec[2]/gigapascal).tolist()
C_micmec_conf6 = (all_C_micmec[6]/gigapascal).tolist()

C_atomic_conf0 = (all_C_atomic[0]/gigapascal).tolist()
C_atomic_conf2 = (all_C_atomic[2]/gigapascal).tolist()
C_atomic_conf6 = (all_C_atomic[6]/gigapascal).tolist()


jdict = {
    "all_K_atomic": [e/gigapascal for e in all_K_atomic],
    "all_E_atomic": [(t[0]/gigapascal, t[1]/gigapascal) for t in all_E_atomic],
    "all_K_micmec": [e/gigapascal for e in all_K_micmec],
    "all_E_micmec": [(t[0]/gigapascal, t[1]/gigapascal) for t in all_E_micmec],
    "K_fcu": K_fcu/gigapascal,
    "K_reo": K_reo/gigapascal,
    "C_micmec_conf0": C_micmec_conf0,
    "C_micmec_conf2": C_micmec_conf2,
    "C_micmec_conf6": C_micmec_conf6,
    "C_atomic_conf0": C_atomic_conf0,
    "C_atomic_conf2": C_atomic_conf2,
    "C_atomic_conf6": C_atomic_conf6,
}

with open("output_configurations.json", "w") as jfile:
    json.dump(jdict, jfile, indent=4)

