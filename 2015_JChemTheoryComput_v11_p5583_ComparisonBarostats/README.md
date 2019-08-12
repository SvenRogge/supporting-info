# Table of contents

This GitHub repository contains the input and raw output data accompanying the manuscript

**A Comparison of Barostats for the Mechanical Characterization of Metal--Organic Frameworks**

by Sven M.J. Rogge, Louis Vanduyfhuys, An Ghysels, Michel Waroquier, Toon Verstraelen, Guillaume Maurin, and Veronique Van Speybroeck. This work was published in *J. Chem. Theory Comput.*, **2015**, *11* (12): 5583--5597 (DOI: 10.1021/acs.jctc.5b00748).


The data presented here is licensed under the CC BY-SA 4.0 international license, a copy of which can be found [here](https://creativecommons.org/licenses/by-sa/4.0/). Under this license, you can copy and redistribute the material in any medium or format as long as you give appropriate credit, provide a link to the license, and indicate if changes were made.

Additional information concerning the data is available upon request from the authors. Please send a mail to Sven.Rogge@UGent.be for more information.

## Software
All simulations were performed with Yaff 1.0, using the module 1.0.barostats11-intel-2014b-Python-2.7.8. The pre- and post-processing Python scripts are written in Python 2.7.




## NPT
The `NPT` folder contains the results and scripts of the constant-temperature constant-pressure simulations performed on MIL-53(Al) at 300 K.

### results
This folder contains the average cell parameters for the three barostats with a barostat relaxation time of 1 ps, for which the trajectories are sliced in a large-pore (lp) phase and a closed-pore (cp) phase. The resulting cell parameters are stored in `cellpars_lp_1ps.csv` and `cellpars_cp_1ps.csv`, respectively. In addition, `max_transitions_all.txt` contains for all barostats, relaxation times, and pressures, the ten simulations in which the lp-to-cp transition occurs the slowest. The tar archive `MIL53_chk.tar.gz` contains initial checkpoint files of structures at a series of volumes, which are used in the `TD_int` folder. In addition, for each of the barostats (subfolders `Berendsen`, `Langevin`, and `NHCMTK`) and each of the barostat time constants a file `transition_Xps.csv` that states at which framenumber an lp-to-cp transition is observed for all simulations that were carried out.

For the subfolder `Berendsen`, the cell parameters of the cp phase at 100 kPa and each of the time constants is stored in `cellpars_Xps_1e5_cp.csv`.

For the subfolder `Langevin`, the average temperature and pressure as well as its standard deviation for the simulations with a 1 ps time constant are stored in `tempress_1ps.csv`.

For the subfolder `NHCMTK`, the average temperature and pressure as well as its standard deviation for the simulations with a 1 ps time constant are stored in `tempress_1ps.csv`. In addition, the frame at which an lp-to-cp transition occurs for each of the time constants is stored in `transition_Xps.csv`. Specifically for a pressure of 1 MPa, additional time constants were employed, and the results are stored in `transition_1MPa_other_masses.csv`. Finally, for the 1 ps relaxation time, the files `transitions_1ps_mXXMPa.csv` contain the frame of transitions for the negative pressures.

### scripts
For each of the barostats (`Berendsen`, `Langevin`, and `NHCMTK`), a subfolder is created that contains a job script `job.sh` as well as an example Yaff input Python script for each of the time constants (`ymd_XXps.csv`). The subfolder `ff` contains an initial checkpoint file (`init.chk`) as well as the QuickFF force field file used here (`pars.txt`).


## TD_int

### results
This folder contains for each of the barostats (subfolders `Berendsen`, `Langevin`, and `NHCMTK`) and each of the barostat time constants a file containing the pressure-versus-volume raw data (`PvsV_Xps_100ps_800Ps.csv`), as well as tar archives containing the bootstrapped free energy profiles  (`bootstrapping_Xps.tar.gz`).

For the subfolder `Berendsen` and the time constants of 5 ps and 10 ps, these files are supplemented by the results obtained for nonconverged simulations (nomenclature as above but ending with `_notconverged`).

For the subfolder `Langevin` and a time constant of 1 ps, the average structures (`av_struct_1ps_XXXXA3.xyz`) and the cell parameters (`cellpars_XX_1ps_XXXXA3.txt`) of the closed-pore (cp), the transition (tr) and large-pore (lp) states are stored.

For the subfolder `NHCMTK`, the file `EPvsV_1ps_100ps_800ps.csv` contains the energy- and pressure-versus-volume data for a time constant of 1 ps. In addition, several bootstrapping tar files have been provided.

### scripts
For each of the barostats (`Berendsen`, `Langevin`, and `NHCMTK`), a subfolder is created that contains a job script `job.sh` as well as an example Yaff input Python script for each of the time constants (`ymd_XXps.csv`). 

