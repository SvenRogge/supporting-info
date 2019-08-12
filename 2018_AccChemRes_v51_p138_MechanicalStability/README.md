# Table of contents

This GitHub repository contains the input and raw output data accompanying the manuscript

**Reliably Modeling the Mechanical Stability of Rigid and Flexible Metal-Organic Frameworks**

by Sven M.J. Rogge, Michel Waroquier, and Veronique Van Speybroeck

This work was published in *Acc. Chem. Res.*, **2018**, *51* (1): 138-148 (DOI: [10.1021/acs.accounts.7b00404](dx.doi.org/10.1021/acs.accounts.7b00404)).

The data presented here is licensed under the CC BY-SA 4.0 international license, a copy of which can be found [here](https://creativecommons.org/licenses/by-sa/4.0/). Under this license, you can copy and redistribute the material in any medium or format as long as you give appropriate credit, provide a link to the license, and indicate if changes were made.

Additional information concerning the data is available upon request from the authors. Please send a mail to Sven.Rogge@UGent.be for more information.

## Software
All molecular dynamics simulations were performed with Yaff 1.0, using the modules 1.0.thesissven2.10.cmm-intel-2015b-Python-2.7.10 (for UiO-66(Zr)) and 1.0.thesissven2.11.cmm-intel-2016a-Python-2.7.11 (for MOF-5(Zn), MIL-47(V), and MIL-53(Al)) if only Yaff was used, except for the $(N,V, \bm \sigma_a =\bm 0, T)$ simulations for MIL-53(Al), which were performed using the module 1.0.barostats12-intel-2014b-Python-2.7.8. 

For those simulations in which Yaff was interfaced with LAMMPS to efficiently calculate the long-range interactions, Yaff 1.0 (1.0.develop.2.14-intel-2015b-Python-2.7.10-HDF5-1.8.15-patch1-serial) and LAMMPS (r12824-intel-2015b) were used. The files `mylammps.py` and `liblammps.py`, provided in the scripts, were needed to achieve this coupling. Parallellization was achieved using mpi4py (module 2.0.0-intel-2015b-Python-2.7.10) and vsc-mympirun (module 3.4.2-intel-2015b-Python-2.7.10-vsc-base-2.4.2). 

The pre- and post-processing Python scripts are written in Python 2.7. 



## UiO66Zr

### results
This folder `UiO66Zr/results` contains the elastic constants at 0 K (`elastic_constants_opt.txt`), the elastic tensors at 300 K and at pressures between 0 MPa and 2000 MPa (`elastic_tensors.tar.gz`), a set of initial structures to construct the pressure-versus-volume equation of state (`structures.tar.gz`), and the resulting pressure-versus-volume equation of state taking into account an equilibration of 50 ps and a production run of 500 ps (`PvsV_eq50ps_av500ps.tar.gz`).

### scripts
The folder `UiO66Zr/scripts/ff` contains the initial structure used for the optimization and constant-pressure simulations (`init.chk`), as well as the QuickFF force field parameter file.

The folder `UiO66Zr/scripts/opt` contains the Python script used to calculate the elastic tensor at 0 K (`StrainStress.py`), as well as the job script used to submit the job (`job.sh`).

The folder `UiO66Zr/scripts/NPT` contains the Python scripts to (re)start the constant-pressure constant-temperature simulations either with (`ymd_yaff_lammps.py` and `ymd_yaff_lammps_restart.py`) or without (`ymd_yaff.py` and `ymd_yaff_restart.py`) LAMMPS. The folder also contains the job scripts used to submit these jobs (`job_xxx.sh`).

The folder `UiO66Zr/scripts/TD_int` contains the Python scripts to (re)start the  $(N,V, \bm \sigma_a =\bm 0, T)$ simulations (`ymd.py` and `ymd_restart.py`) as well as the job scripts used to submit these jobs (`job.sh` and `job_restart.sh`).

The folder `UiO66Zr/scripts/LAMMPS` contains the Python scripts to generate the LAMMPS table containing the long-range interactions (`gentable.py`), the resulting table (`lammps_smoothei2.table`), and the two files needed to interface Yaff and LAMMPS (`liblammps.py` and `mylammps.py`) for the specified modules of these software codes.


## MOF5Zn

### results
This folder `MOF5Zn/results` contains the elastic constants at 0 K (`elastic_constants_opt.txt`), the elastic tensors at 300 K and at pressures between -500 MPa and 180 MPa (`elastic_tensors.tar.gz`), a set of initial structures to construct the pressure-versus-volume equation of state (`structures.tar.gz`), and the resulting pressure-versus-volume equation of state taking into account an equilibration of 100 ps and a production run of 500 ps (`PvsV_eq100ps_av500ps.tar.gz`).

### scripts
The folder `MOF5Zn/scripts/ff` contains the initial structure used for the optimization and constant-pressure simulations (`init.chk`), as well as the QuickFF force field parameter file.

The folder `MOF5Zn/scripts/opt` contains the Python script used to calculate the elastic tensor at 0 K (`StrainStress.py`), as well as the job script used to submit the job (`job.sh`).

The folder `MOF5Zn/scripts/NPT` contains the Python scripts to (re)start the constant-pressure constant-temperature simulations (`ymd.py` and `ymd_restart.py`) as well as the job scripts used to submit these jobs (`job.sh` and `job_restart.sh`).

The folder `MOF5Zn/scripts/TD_int` contains the Python scripts to (re)start the  $(N,V, \bm \sigma_a =\bm 0, T)$ simulations (`ymd.py` and `ymd_restart.py`) as well as the job scripts used to submit these jobs (`job.sh` and `job_restart.sh`).

The folder `MOF5Zn/scripts/LAMMPS` contains the Python scripts to generate the LAMMPS table containing the long-range interactions (`gentable.py`), the resulting table (`lammps_smoothei2.table`), and the two files needed to interface Yaff and LAMMPS (`liblammps.py` and `mylammps.py`) for the specified modules of these software codes.


## MIL47V_lp

### results
This folder `MIL47V_lp/results` contains the elastic constants at 0 K (`elastic_constants_opt.txt`), the elastic tensors at 300 K and at pressures between 0 MPa and 105 MPa (`elastic_tensors.tar.gz`), and the resulting pressure-versus-volume equation of state taking into account an equilibration of 100 ps and a production run of 500 ps (`PvsV_eq100ps_av500ps.tar.gz`).

### scripts
The folder `MIL47V_lp/scripts/ff` contains the initial structure used for the optimization and constant-pressure simulations (`init.chk`), as well as the QuickFF force field parameter file.

The folder `MIL47V_lp/scripts/opt` contains the Python script used to calculate the elastic tensor at 0 K (`StrainStress.py`), as well as the job script used to submit the job (`job.sh`).

The folder `MIL47V_lp/scripts/NPT` contains the Python scripts to (re)start the constant-pressure constant-temperature simulations (`ymd.py` and `ymd_restart.py`) as well as the job scripts used to submit these jobs (`job.sh` and `job_restart.sh`).

The folder `MIL47V_lp/scripts/LAMMPS` contains the Python scripts to generate the LAMMPS table containing the long-range interactions (`gentable.py`), the resulting table (`lammps_smoothei2.table`), and the two files needed to interface Yaff and LAMMPS (`liblammps.py` and `mylammps.py`) for the specified modules of these software codes.


## MIL53Al_lp

### results
This folder `MIL53Al_lp/results` contains the elastic constants at 0 K (`elastic_constants_opt.txt`), the elastic tensors at 300 K and at pressures between -200 MPa and -10 MPa (`elastic_tensors.tar.gz`), a set of initial structures to construct the pressure-versus-volume equation of state (`structures.tar.gz`), and the resulting pressure-versus-volume equation of state taking into account an equilibration of 100 ps and a production run of 700 ps (`PvsV_eq100ps_av700ps.tar.gz`).

### scripts
The folder `MIL53Al_lp/scripts/ff` contains the initial structure used for the optimization and constant-pressure simulations (`init.chk`), as well as the QuickFF force field parameter file.

The folder `MIL53Al_lp/scripts/opt` contains the Python script used to calculate the elastic tensor at 0 K (`StrainStress.py`), as well as the job script used to submit the job (`job.sh`).

The folder `MIL53Al_lp/scripts/NPT` contains the Python scripts to (re)start the constant-pressure constant-temperature simulations (`ymd.py` and `ymd_restart.py`) as well as the job scripts used to submit these jobs (`job.sh` and `job_restart.sh`).

The folder `MIL53Al_lp/scripts/TD_int` contains the Python scripts to start the  $(N,V, \bm \sigma_a =\bm 0, T)$ simulations (`ymd.py`) as well as the job script used to submit these jobs (`job.sh`).

The folder `MIL53Al_lp/scripts/LAMMPS` contains the Python scripts to generate the LAMMPS table containing the long-range interactions (`gentable.py`), the resulting table (`lammps_smoothei2.table`), and the two files needed to interface Yaff and LAMMPS (`liblammps.py` and `mylammps.py`) for the specified modules of these software codes.


## MIL53Al_cp

### results
This folder `MIL53Al_cp/results` contains the elastic constants at 0 K (`elastic_constants_opt.txt`) and the elastic tensors at 300 K and at pressures between -100 MPa and 100 MPa (`elastic_tensors.tar.gz`).

### scripts
The folder `MIL53Al_cp/scripts/ff` contains the initial structure used for the optimization and constant-pressure simulations (`init.chk`), as well as the QuickFF force field parameter file.

The folder `MIL53Al_cp/scripts/opt` contains the Python script used to calculate the elastic tensor at 0 K (`StrainStress.py`), as well as the job script used to submit the job (`job.sh`).

The folder `MIL53Al_cp/scripts/NPT` contains the Python scripts to (re)start the constant-pressure constant-temperature simulations (`ymd.py` and `ymd_restart.py`) as well as the job scripts used to submit these jobs (`job.sh` and `job_restart.sh`).

The folder `MIL53Al_cp/scripts/LAMMPS` contains the Python scripts to generate the LAMMPS table containing the long-range interactions (`gentable.py`), the resulting table (`lammps_smoothei2.table`), and the two files needed to interface (`liblammps.py` and `mylammps.py`) for the specified modules of these software codes.
