# Table of contents

This GitHub repository contains the input and raw output data accompanying the manuscript

**The Importance of Cell Shape Sampling To Accurately Predict Flexibility in Metal--Organic Frameworks**

by Sven M.J. Rogge, Senne Caroes, Ruben Demuynck, Michel Waroquier, Veronique Van Speybroeck, and An Ghysels.

This work was published in *J. Chem. Theory Comput.*, **2018**, *14* (3): 1186-1197 (DOI: [10.1021/acs.jctc.7b01134](https:/dx.doi.org/10.1021/acs.jctc.7b01134)).

The data presented here is licensed under the CC BY-SA 4.0 international license, a copy of which can be found [here](https://creativecommons.org/licenses/by-sa/4.0/). Under this license, you can copy and redistribute the material in any medium or format as long as you give appropriate credit, provide a link to the license, and indicate if changes were made.

Additional information concerning the data is available upon request from the authors. Please send a mail to Sven.Rogge@UGent.be for more information.


# Software
All molecular dynamics simulations were performed with Yaff 1.0, using the module 1.0.develop.2.15-intel-2015b-Python-2.7.10-HDF5-1.8.15-patch1-serial. 

For the MIL-53(Al) simulations, Yaff was interfaced with the LAMMPS module r12824-intel-2015b. The files `mylammps.py` and `liblammps.py`, provided in the scripts, were needed to achieve this coupling. Parallelization was achieved using mpi4py (module 2.0.0-intel-2015b-Python-2.7.10) and vsc-mympirun (module 3.4.2-intel-2015b-Python-2.7.10-vsc-base-2.4.2). 

For the DUT-49(Cu) simulations, Yaff was interfaced with DL_POLY version 2.20.

The pre- and post-processing Python scripts are written in Python 2.7. 


# MIL53Al

## results
This folder contains the necessary information extracted from the MD trajectories to create the pressure and free energy equations of state using the `MTD`, `TI_full`, `TI_ipol1`, `TI_ipol2`, `TI_mean`, and `TI_snap` methods above in `profile.tar.gz` (Figure 3(a), Figure 4(a), Figure 5, and Figure 6 in the article). For the `TI_ipol1` and `TI_ipol2` methods, it also contains the generalized pressure equation of state, reported in Figure 4(b) in the article. For the `TI_full` method, the results obtained with the 1x2x1 and 2x4x2 supercells are reported separately. These data, in `.npy` format, can be visualized using the scripts provided in `scripts/analysis`.

## scripts
The folder `analysis` contains all post-analysis scripts:

 - `2D_hist.py` plots the histogram between the volume and the components of the cell shape (see Figure S2)
 - `average_cell.py`determines from an $(N, V, \bm \sigma_a = \bm 0, T)$ simulation the average cell shape, which is then used for the $(N, V, \mathbf{h}_0, T)$ simulations in the `TI_mean` method
 - `construct_profile_MTD.py` and `construct_profile_TI_full.py` are two example scripts used to convert the MD trajectories to the `.npy` data for the simulations based on metadynamics and thermodynamic integration, respectively
 - `plot_data_npy.py` reads in all data from `results` to generate the figures shown in the article

The folder `ff` contains the Yaff force field file (`pars.txt`), the script necessary to generate the LAMMPS table of non-covalent interactions (`gentable_pars.py`), and the resulting LAMMPS table (`lammps_smoothei2.table`).

The folder `LAMMPS` contains the files `liblammps.py` and `mylammps.py` necessary for the coupling between Yaff and LAMMPS.

The folders `MTD`, `TI_full`, `TI_ipol1`, `TI_ipol2`, and `TI_snap` contain all necessary input files for these simulations, including the `ymd.py` files to start the Yaff simulation and the corresponding `job.sh` job script. For the `MTD` folder, also the definition of the hills (`colvar_senne.py`), the MTD module (`DMTD_senne.py` and `DMTD_senne_restart.py`), and the initial `.chk` files for the ten independent simulations (`init_random_chks.tar.gz`) are provided. The job scripts for all thermodynamic integration simulations are similar; hence only the ones for `TI_full` are provided. For all thermodynamic integration methods except for `TI_full` the initial series of `.chk` files is provided in `chk_files.tar.gz`.


# DUT49Cu

## results
The folder `DUT49Cu/results` contains the force-field optimized cp state in the subfolder `relaxation` as well as the necessary information extracted from the MD trajectories to create the pressure and free energy equations of state using the `TI_full` and `TI_ipol` methods above (Figure 3(b) and Figure 7 in the article) in `profile.tar.gz`.  These data, in `.npy` format, can be visualized using the scripts provided in `scripts/analysis`.

## scripts
The folder `analysis` contains the post-analysis script `plot_data_npy.py`, which reads in all data from `results` to generate the figures shown in the article.

The folder `ff` contains the Yaff and DL_POLY force field and input files, as well as python scripts to convert between the two.

The folder `LAMMPS` contains the files `liblammps.py` and `mylammps.py` necessary for the coupling between Yaff and LAMMPS.

The folders `NPT`, `relaxation`, `TI_full` and `TI_ipol` contain the `ymd.py` files to start the Yaff simulation with the DL_POLY force field.
