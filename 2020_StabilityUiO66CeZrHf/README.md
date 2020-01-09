# Table of contents

This GitHub repository contains the input and raw output data accompanying the manuscript

**Systematically Isolating the Intrinsic and Defect-Controlled High-Pressure Mechanical Stability of Bimetallic UiO-66 Materials**

by Sven M. J. Rogge, Pascal G. Yot, Jannick Jacobsen, Francesco Muniz-Miranda, Steven Vandenbrande, Jonas Gosch, Vanessa Ortiz, Inès E. Collings, Sabine Devautour-Vinot, Guillame Maurin, Norbert Stock, and Veronique Van Speybroeck.

The data presented here is licensed under the CC BY-SA 4.0 international license, a copy of which can be found [here](https://creativecommons.org/licenses/by-sa/4.0/). Under this license, you can copy and redistribute the material in any medium or format as long as you give appropriate credit, provide a link to the license, and indicate if changes were made.

Additional information concerning the data is available upon request from the authors. Please send a mail to Sven.Rogge@UGent.be for more information.


## Software
*Gaussian*
The 22 isolated inorganic building blocks have been relaxed using the Becker 3-parameters hybrid exchange functional in combination with the Lee-Yang-Parr correlation (B3LYP) functional. The split-valence douyble-zeta 6-31G(d,p) basis set has been employed to describe the orbitals of the hydrogen, carbon, and oxygen atoms, whereas the Stuttgart combined basis set and pseudopotential (the so-called 'SDD' basis set) has been used to model the metal atoms (cerium, hafnium, and zirconium). All these methods have been adopted as implemented in the Gaussian 16 software suite, using the module g16_E.01-intel-2017a. The standard convergence criteria were adopted for the optimization of the geometries (corresponding to forces and atomic displacements below the 2.72 x 10$^{-4}$ eV/bohr and 4 x 10$^{-5}$ bohr thresholds, respectively).

*QuickFF and Horton*
To generate force fields for each of the 22 periodic UiO-66 structures, desktop-versions of QuickFF v2 and HORTON were used.

*Yaff and LAMMPS*
All molecular dynamics simulations were performed with Yaff 1.4.2, using the module 1.4.2-intel-2017b-Python-2.7.14. For those simulations in which Yaff was interfaced with LAMMPS to efficiently calculate the long-range interactions, the r12824-intel-2017b LAMMPS module was used. The files `mylammps.py` and `liblammps.py`, provided in the scripts, were needed to achieve this coupling. Parallellization was achieved using mpi4py (module 2.0.0-intel-2017b-Python-2.7.14).

*Pre- and post-processing scripts*
The pre- and post-processing Python scripts are written in Python 2.7. To convert the trajectory files to crystallographic information files (CIFs), CON3F was used (module version_CIF2-intel-2015a-Python-2.7.10).

## construct_ff_quickff_horton

### results
The folder `construct_ff_quickff_horton/results` contains the force-field optimized structure files for each of the 22 periodic materials (`init.tar.gz`) and the corresponding periodic force fields (`pars.tar.gz`)


## NPT_yaff

### results
This folder contains 22 sets of initial structures to construct the pressure-versus-volume equations of state for each of the 22 materials in the paper (`structures_xxx.tar.gz`).

### scripts
The file `extract_structures.py` extracts from a given constant-pressure constant-temperature simulation the structures closest to a predefined volume grid and within a given error threshold.

The `job.sh` file is an example of the files used to submit the molecular dynamics simulations to the HPC.

The `ymd.py` file  is an example of the input scripts used to start the molecular dynamics simulations.

The file `ff/ff_pars.tar.gz` contains the force field parameter files for each of the 22 materials.

The file `ff/initial_structures.tar.gz` contains the initial structures used to start the constant-temperature constant-pressure simulations.


## TD_int_yaff_lammps

### results
This folder contains the crystallographic information files (CIFs) for the averaged structures at all volumes of all 22 materials (`CIFs.tar.gz`) and the pressure-versus-volume equations of state for all materials taking into account an equilibration of 100 ps and a production run of 900 ps (`PvsV_eq100ps_av900ps.tar.gz`).

### scripts

The `job.sh` file is an example of the files used to submit the molecular dynamics simulations to the HPC.

The files `liblammps.py` and `mylammps.py` are needed to interface Yaff and LAMMPS.

The file `PvsV.py` was used to construct the pressure-versus-volume equations of state, taking an equilibration of 50 ps into account.

The `ymd.py` file  is an example of the input scripts used to start the molecular dynamics simulations.

The file `ff/ff_pars_yaff.tar.gz` contains the Yaff force field parameter files for each of the 22 materials.

The file `ff/gentable.py`was used to construct the LAMMPS tables to efficiently calculate the long-range interactions.

The file `ff/lammps_smoothei2.tar.gz` contains the tabulated values for the long-range interactions used by LAMMPS for each of the 22 materials.

The files `ff/unique_ffatypes_2_Ce.py` and `ff/unique_ffatypes_2_Hf.py` were used to automatically redefine the force field atom types to distinguish between oxygen and carbon atoms that have a different charge.

The file `ff/ymd_test.py`was used to test whether the Yaff and LAMMPS energies are in agreement.
