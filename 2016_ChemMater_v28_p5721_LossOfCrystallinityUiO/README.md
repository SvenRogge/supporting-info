# Table of contents

This GitHub repository contains the input and raw output data accompanying the manuscript

**Thermodynamic Insight in the High-Pressure Behavior of UiO-66: Effect of Linker Defects and Linker Expansion**

by Sven M.J. Rogge, Jelle Wieme, Louis Vanduyfhuys, Steven Vandenbrande, Guillaume Maurin, Toon Verstraelen, Michel Waroquier, and Veronique Van Speybroeck

This work was published in *Chem. Mater.*, **2016**, *28* (16): 5721-5732 (DOI: [10.1021/acs.chemmater.6b01956](https://dx.doi.org/10.1021/acs.chemmater.6b01956)).

The data presented here is licensed under the CC BY-SA 4.0 international license, a copy of which can be found [here](https://creativecommons.org/licenses/by-sa/4.0/). Under this license, you can copy and redistribute the material in any medium or format as long as you give appropriate credit, provide a link to the license, and indicate if changes were made.

Additional information concerning the data is available upon request from the authors. Please send a mail to Sven.Rogge@UGent.be for more information.

## Software
All molecular dynamics simulations were performed with Yaff 1.0, using the module 1.0.thesissven2.10.cmm-intel-2015b-Python-2.7.10 if only Yaff was used. For those simulations in which Yaff was interfaced with LAMMPS to efficiently calculate the long-range interactions, Yaff 1.0 (module 1.0.develop.2.14-intel-2015b-Python-2.7.10-HDF5-1.8.15-patch1-serial) and LAMMPS (module r12824-intel-2015b) were used. The files `mylammps.py` and `liblammps.py`, provided in the scripts, were needed to achieve this coupling. Parallellization was achieved using mpi4py (module 2.0.0-intel-2015b-Python-2.7.10) and vsc-mympirun (module 3.4.2-intel-2015b-Python-2.7.10-vsc-base-2.4.2). 

The pre- and post-processing Python scripts are written in Python 2.7. To convert the trajectory files to crystallographic information files (CIFs), CON3F was used (module version_CIF2-intel-2015a-Python-2.7.10).



## NPT_yaff_lammps
The folder `NPT_yaff_lammps` contains the data obtained by constant-temperature constant-pressure simulations for all UiO-type materials.

### results
This folder contains the elastic tensors of the defect-free UiO-66 structure at pressures between 0 MPa and 2000 MPa (`elastic_tensors.tar.gz`) as well as a set of initial structures to construct the pressure-versus-volume equations of state for each of the eleven materials in the paper (`structures_xxx.tar.gz`).

### scripts
The different `job_xxx.sh` files were used to submit the molecular dynamics simulations to the HPC, either with Yaff alone or with the Yaff+LAMMPS interface.

The file `lammps_smoothei2.table` contains the tabulated values for the long-range interactions used by LAMMPS.

The files `liblammps.py` and `mylammps.py` are needed to interface Yaff and LAMMPS.

The different `ymd_xxx.py` files contain the input scripts to start the molecular dynamics simulations, either with Yaff alone or with the Yaff+LAMMPS interface.

The folder `ff` contains the force field parameter files as well as the the initial structures used to start the constant-temperature constant-pressure simulations for each of the eleven materials .


## TD_int_yaff
The folder `TD_int_yaff_lammps` contains all data obtained by constant-temperature constant-volume simulations in the $(N,V, \bm \sigma_a =\bm 0, T)$ ensemble for all UiO-type materials.

### results
This folder contains the crystallographic information files (CIFs) for the averaged structures of all defect-free materials (`CIFs_UiO6x.tar.gz`), the equilibrium structures for each of the eleven materials (`equilibrium_structures_tar.gz`), and the pressure-versus-volume equations of state for all materials taking into account an equilibration of 50 ps and a production run of 500 ps (`PvsV_eq50ps_av500ps.tar.gz`).

### scripts

The different `job_xxx.sh` files were used to submit the molecular dynamics simulations to the HPC.

The file `make_cif.sh` was used to convert the trajectory files to a CIF of the averaged structure using CON3F.

The file `PvsV.py` was used to construct the pressure-versus-volume equations of state, taking an equilibration of 50 ps into account.

The different `ymd_xxx.py` files contain the input scripts to start the molecular dynamics simulations.
