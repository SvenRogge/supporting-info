# Table of contents

This GitHub repository contains the input and raw output data accompanying the manuscript

**"Structure-Mechanical Stability Relations of Metal-Organic Frameworks via Machine Learning"**

by Peyman Z. Moghadam, Sven M.J. Rogge, Aurelia Li, Chun-Man Chow, Jelle Wieme, Noushin Moharrami, Marta Aragones-Anglada, Gareth Conduit, Diego A. Gomez-Gualdron, Veronique Van Speybroeck, and David Fairen-Jimenez.

This work was published in *Matter*, **2019**, *1* (1): 219--234 (DOI: [10.1016/j.matt.2019.03.002](https://www.cell.com/matter/fulltext/S2590-2385(19)30006-2)).

The data presented here is licensed under the CC BY-SA 4.0 international license, a copy of which can be found [here](https://creativecommons.org/licenses/by-sa/4.0/). Under this license, you can copy and redistribute the material in any medium or format as long as you give appropriate credit, provide a link to the license, and indicate if changes were made.

Additional information concerning the data is available upon request from the authors. Please send a mail to Sven.Rogge@UGent.be for more information.


## Software
*Gaussian*
The isolated metal brick and fourteen isolated ligands were optimized with the Gaussian16 software suite (g16_E.01-intel-2017a) using the Becker 3-parameters hybrid exchange functional in combination with the Lee-Yang-Parr correlation (B3LYP) functional. The 6-311G(d,p) Pople basis set was used for all atoms, except for zirconium, which was described using the LANL2DZ basis set, including an effective core potential. The standard convergence criteria were adopted for the optimization of the geometries (corresponding to forces and atomic displacements below the 2.72 x 10<sup>-4</sup> eV/bohr and 4 x 10<sup>-5</sup> bohr thresholds, respectively). A vibrational frequency analysis ensured that the optimized structures correspond to minima on the potential energy surface.

*GPAW and HORTON*
The GPAW module 0.11.0.13004-intel-2015b-Python-2.7.10 and the corresponding GPAW-setups module 0.9.11271-linux-x86_64 have been employed to determine the charges using the PBE exchange-correlation functional and a grid spacing of 0.15 &#8491;. The resulting charges were then extracted using the MBIS scheme as implemented in HORTON. 

*QuickFF*
From the Hessian and the MBIS charges of the cluster models determined above, flexible and _ab initio_-based force fields for the isolated cluster models were derived. The covalent part of the force field, which contains diagonal terms that describe bonds, bends, out-of-plane distances, and torsion angles was derived using the in-house QuickFF software package (version 2.2.0). The van der Waals interactions in these force field models were added *a posteriori* based on the MM3 parameters. These force field models of the isolated building blocks were combined to obtain force field models for the materials discussed in the manuscript.

*Yaff and LAMMPS*
All force field optimizations and molecular dynamics simulations were performed with Yaff. For the optimizations and the subsequent extraction of the elastic tensor at 0 K, the module 1.4.2-intel-2017b-Python-2.7.14 was used. For the (*N*, *P*, **&sigma;**<sub>a</sub> = **0**, *T*) and  (*N*, *V*, **&sigma;**<sub>a</sub> = **0**, *T*)  simulations at 300 K, the module 1.1.3-intel-2017a-Python-2.7.13 was used and, if applicable, interfaced with LAMMPS (module r12824-intel-2017a) to efficiently calculate the long-range interactions. The files `mylammps.py` and `liblammps.py`, provided in the scripts, were needed to achieve this coupling. Parallellization was achieved using mpi4py (module 2.0.0-intel-2017a-Python-2.7.13-timed-pingpong) and vsc-mympirun (module 4.0.2). 

*Pre- and post-processing scripts*
The pre- and post-processing Python scripts are written in Python 2.7.


## construct_ff_quickff_horton_gaussian

### input
This folder contains the initial geometries for the inorganic brick and fourteen organic ligands which were used to initialize the Gaussian optimizations (`init_brick_X.com`). It also contains the QuickFF input file (`quickffrc`) with the parameters used to derive the force fields.

### results
This folder contains two subfolders.  The **gaussian** subfolder contains the optimized structures (`opt_*.xyz`) for the isolated building blocks. The **quickff** subfolder contains the force fields (both the covalent part,`pars_*_cov.txt`, and the electrostatic part, `pars_*_ei.txt`) and Yaff system files (`init_*.chk`) for these isolated building blocks, which were combined to yield the force fields for the periodic structures discussed in the manuscript.

### scripts
This folder contains the job scripts to start the Gaussian optimization (`gaussian_opt.sh`), the Gaussian frequency analysis (`gaussian_freq.sh`), the GPAW calculation (`gpaw.sh`), and the extraction of the MBIS charges (`run_denspart.sh`). It also contains the file `make_sample.py`, which is used to add the correct atom types and bonds to the different building blocks.


## opt_yaff

### input
This folder contains the force fields (`pars_*.txt`) and Yaff system files (`init_*.chk`) for the fourteen periodic materials (UiO-L1 to UiO-L14) discussed in the manuscript. These initial structures were used to start the optimizations.

### results
This folder contains the 0 K elastic constants for the fourteen materials (UiO-L1 to UiO-L14) discussed in the manuscript, as determined using the derived force fields (`etens_0K_*.txt`). 

### scripts
This folder contains the file `elastic_constants.py`, which optimizes a given geometry and afterwards derives the elastic constans in equilibrium. It also contains the corresponding job script (`elastic_constants.sh`).



## NPT_yaff_lammps

### input
This folder contains the force fields (`pars_*.txt`) and Yaff system files (`init_*.chk`) for the fourteen materials discussed in the manuscript (UiO-L1 to UiO-L14). These initial structures were used to start the (*N*, *P*, **&sigma;**<sub>a</sub> = **0**, *T*)  simulations. 

### results
This folder contains the elastic tensors of the fourteen materials discussed in the manuscript  (UiO-L1 to UiO-L14) at various pressures below the instability pressure (`elastic_tensors_UiO_*.tar.gz`). It also contains a set of initial structures to construct the pressure-versus-volume equations of state for each of these fourteen materials discussed in the paper (`structures_UiO-*.tar.gz`).

### scripts
The file `add_cnt.py` was needed to be run after completion of each h5-file to allow for the restart of this simulation.

The file `calc_elast.py` calculates the elastic constants for a cubic material.

The file `extract_structures.py` extracts from a given constant-pressure constant-temperature simulation the structures closest to a predefined voume grid and within a given error threshold.

The file `gentable.py` generates the LAMMPS table for the long-range interactions starting from the Yaff force field and system object.

The different `job_xxx.sh` files were used to submit the molecular dynamics simulations, either with Yaff alone or with the Yaff+LAMMPS interface (original submission and restart).

The files `liblammps.py` and `mylammps.py` are needed to interface Yaff and LAMMPS.

The job scripts `makedirs_*.sh` and `restartdirs_yaff_lammps.sh` were used to prepare the various directories.

The different `ymd_xxx.py` files contain the input scripts to start the molecular dynamics simulations, either with Yaff alone or with the Yaff+LAMMPS interface.


## TD_int_yaff

### results
This folder contains the pressure-versus-volume equations of state for the fourteen materials discussed in the manuscript (UiO-L1 to UiO-L14), taking into account an equilibration of 100 ps and a production run of 900 ps (`PvsV_UiO-*_eq100ps_av900ps.csv`).

### scripts
The file `add_cnt.py` was needed to be run after completion of each h5-file to allow for the restart of this simulation.

The different `job_xxx.sh` files were used to submit the molecular dynamics simulations (original submission and restart).

The file `PvsV.py` was used to construct the pressure-versus-volume equations of state, taking an equilibration of 100 ps into account.

The different `ymd_xxx.py` files contain the input scripts to start the molecular dynamics simulations.
