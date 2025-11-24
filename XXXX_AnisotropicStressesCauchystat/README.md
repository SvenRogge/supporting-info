


# Table of contents

This GitHub repository contains the input data accompanying the manuscript

**Challenges in modelling anisotropic stresses in soft polymorphic materials**

by Jelto Neirynck, Sander Geerinckx, and Sven M. J. Rogge.

The data presented here is licensed under the CC BY-SA 4.0 international license, a copy of which can be found [here](https://creativecommons.org/licenses/by-sa/4.0/). Under this license, you can copy and redistribute the material in any medium or format as long as you give appropriate credit, provide a link to the license, and indicate if changes were made.

Additional information concerning the data is available upon request from the authors. Please email Sven.Rogge@UGent.be for more information.

## Software
All simulations were performed with LAMMPS (feature release 4 February 2025), which is available for free download at https://www.lammps.org/download.html. The conversion of the YAFF input files and force fields to LAMMPS format used YAFF version 1.6.0.


## Archive

This archive contains two folders: *structures* and *input_files*.

### Structures
This folder contains `xyz` files for both the MIL-53(Al) (3x6x3 supercell) and COF-5 (4x4x4 supercell) materials, for easy access to the structures.

### Input_files
This folder contains the actual input scripts used in this manuscript. The input files are grouped by structure into two main folders: *MIL-53* and *COF-5*.

Each final subfolder in *MIL-53* and *COF-5*, excluding the *yaff_input* folder, contains four input files:

* `lammps.in`: contains all MD settings for the simulation run. To adjust the magnitude of the stress, one should update the number **200** in `lammps.in`, *e.g.*, in the line `fix * all npt/cauchy temp 300.0 300.0 100.0 x 200 200 1000 y 0 0 1000 z 0 0 1000 xy 0 0 1000 xz 0 0 1000 yz 0 0 1000 alpha 0.01 nreset 10` to the desired pressure in atm. The parameters *alpha* and *nreset* can be altered by changing the values **0.01** and **10**, respectively. The seed can be adjusted by changing the number **1** in the line `velocity all create 300.0 1`. 
* `lammps.table`: contains the structure and all the bonded interactions
* `eibonded.table`: contains all electrostatic interactions
* `lammps.table`: contains the van der Waals part of the forcefield. 

To start a simulation, you can call the LAMMPS executable on `lammps.in` in the desired folder.

The subfolders are as follows:

#### MIL-53

The *MIL-53* folder contains three subfolders. 

1. *Cauchystat* contains all input scripts for simulations performed with the Cauchystat in the four directions: *x-tensile*, *z-tensile*, *x-compressive*, and *x-tensile*.
2. *MTTK* contains all input scripts using the Martyna-Tuckerman-Tobias-Klein (MTTK) barostat to simulate hydrostatic conditions for the selection of the proper MIL-53 simulation cell size: *1x2x1*, *2x4x2*, *3x6x3*, and *4x8x4*.
3. *yaff_input* contains the original starting structure and QuickFF force field in the format used by YAFF. In this manuscript, they were converted to LAMMPS.

In addition, it contains the file `yaff_to_lammps.py` used to convert the YAFF input structure and force field to LAMMPS.

#### COF-5
*COF-5* contains, next to the *yaff_input* folder discussed above for MIL-53(Al), three subfolders. They correspond with the three directions in which the shear stress simulations were performed: *xz-shear*, *yz-shear*, and *30_degree-shear*.

In addition, it contains the file `yaff_to_lammps.py` used to convert the YAFF input structure and force field to LAMMPS.
