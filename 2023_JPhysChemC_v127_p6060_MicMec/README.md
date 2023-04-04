# Table of contents

This GitHub repository contains the input and raw output data accompanying the manuscript

**MicMec: Developing the Micromechanical Model to Investigate the Mechanics of Correlated Node Defects in UiO-66**

by Joachim Vandewalle, Juul S. De Vos, and Sven M.J. Rogge.

This work was published in *J. Phys. Chem. C*, **2023**, *127* (12): 6060-6070 (DOI: [10.1021/acs.jpcc.3c00451](https://dx.doi.org/10.1021/acs.jpcc.3c00451)).

The data presented here is licensed under the CC BY-SA 4.0 international license, a copy of which can be found [here](https://creativecommons.org/licenses/by-sa/4.0/). Under this license, you can copy and redistribute the material in any medium or format as long as you give appropriate credit, provide a link to the license, and indicate if changes were made.

Additional information concerning the data is available upon request from the authors. Please send a mail to Sven.Rogge@UGent.be for more information.

## Software

The MicMec software is freely available at https://github.com/Jlvdwall/micmec, with
documentation at https://jlvdwall.github.io/micmec/. 


## Archive

This archive consists of two main folders. The folder **test_system** contains all data pertaining to the micromechanical simulations on the test system discussed in sections 3.1 and 3.2 of the manuscript. The folder **UiO-66** contains all data pertaining to the micromechanical simulations on the **fcu**, **reo**, and **fcu**:**reo** models discussed in section 3.3 of the manuscript. 

The structure of each of these folders is discussed below.

### test_system/node_deformation

This subfolder describes the data underlying Figure 4 of the manuscript. Specifically, the **input** subfolder contains:
- `2x2x2_test_micmec.chk` : 
    A 2 x 2 x 2, micromechanical, periodic test system, complete with coordinates and coarse-grained parameters.
- `type_test_micmec.pickle` : 
    The pickle file containing the properties of the test type, which is used in the Micromechanical Model Builder to help create the test system.
- `static_scan.py` : 
    Python script that performs a static scan of a micromechanical system. Use with command-line arguments (see `--help` for more information).
    
In turn, the **output** subfolder contains:
- `output_energy.json` and `fig_energy.py` :
The output data underlying Figure 4b, which can be visualized using the python script.
- `output_forces.json` and `fig_forces.py` :
The output data underlying Figure 4c, which can be visualized using the python script.

### test_system/volumetric_compression

This subfolder  describes the data underlying Figure 5 of the manuscript and Figure S1 of the Supporting Information. Specifically, the **input** subfolder contains:
- `2x2x2_test_micmec.chk` : 
    A 2 x 2 x 2, micromechanical, periodic test system, complete with coordinates and coarse-grained parameters.
- `3x3x3_test_micmec.chk` : 
    A 3 x 3 x 3, micromechanical, periodic test system, complete with coordinates and coarse-grained parameters.
- `type_test_micmec.pickle` : 
    The pickle file containing the properties of the test type, which is used in the Micromechanical Model Builder to help create the test system.
- `relaxed_scan.py` : 
    Python script that performs a relaxed scan of a micromechanical system. Use with command-line arguments (see `--help` for more information).
    
In turn, the **output** subfolder contains:
- `output_relaxed_scan_new.json` :
The output data underlying Figure 5 and Figure S1b for both the 2 x 2 x 2 and 3 x 3 x 3 models, obtained using Eq. S1 for the elastic deformation energy.
- `output_relaxed_scan_old.json` :
The output data underlying Figure S1afor both the 2 x 2 x 2 and 3 x 3 x 3 models, obtained using Eq. S2 for the elastic deformation energy.
- `fig_relaxed_scan.py` :
Python script used to visualize the two aforementioned data files.

### test_system/timestep
This subfolder describes the data underlying Figure 6 of the manuscript. Specifically, the **input** subfolder contains:
- `2x2x2_test_micmec.chk` : 
    A 2 x 2 x 2, micromechanical, periodic test system, complete with coordinates and coarse-grained parameters.
- `type_test_micmec.pickle` : 
    The pickle file containing the properties of the test type, which is used in the Micromechanical Model Builder to help create the test system.
- `md_timesteps.py` : 
    Python script that performs a set of micromechanical simulations, differing in the timestep used. Use with command-line arguments (see `--help` for more information).

In turn, the **output** subfolder contains:
- `timestepX.h5` :
The output data underlying Figure 6 for each of the five different timesteps `X` reported there..
- `fig_relaxed_scan.py` :
Python script used to visualize the five aforementioned data files.

### UiO-66/input
This subfolder contains the input data for all UiO-66 simulations. Specifically, the **atomic** subfolder contains:
- `3x3x3_confX_atomic.chk` :
 The 3 x 3 x 3, atomic system for each of the eight configurations `X` visualized in Figure 7, as well as for the pure **fcu** (`X=0`) and the pure **reo** (`X=9`) structure
 - `fcu.chk` and `reo.chk` :
 The 1 x 1 x 1, atomic system for the pure **fcu** and pure **reo** structures
 - `pars_fcu_atomic.txt` and `pars_reo_atomic.txt` :
 The force field parameters for the **fcu** and **reo** structures, used for all aforementioned systems.
 - `stress_strain_atomic.py` :
 Python script to determine the elastic stiffness tensors of the atomic systems by imposing a strain and measuring the resulting strain.
 - `supercell.py` :
 Python script used to generate the `3x3x3_confX_atomic.chk` files for the different **fcu**:**reo** configurations. 
 
In turn, the **micromechanical** subfolder contains:
- `3x3x3_confX_micmec.chk` :
 The 3 x 3 x 3, micromechanical system for each of the eight configurations `X` visualized in Figure 7, as well as for the pure **fcu** (`X=0`) and the pure **reo** (`X=9`) structure
 - `type_1x1x1_fcu_atomic.pickle` and `type_1x1x1_reo_atomic.pickle` :
 The elastic constants obtained from the 1 x 1 x 1, atomic **fcu** and **reo** structures, which are used as force field parameters in the micromechanical simulations.
 - `stress_strain_micromechanical.py` :
 Python script to determine the elastic stiffness tensors of the micromechanical systems by imposing a strain and measuring the resulting strain.
 
### UiO-66/output
 
 This subfolder contains the output data for all UiO-66 simulations:
 - `type_3x3x3_confX_atomic.pickle` :
Basic output information, including the elastic stifness tensor, of the 3 x 3 x 3, atomic system for each of the eight configurations `X` visualized in Figure 7, as well as for the pure **fcu** (`X=0`) and the pure **reo** (`X=9`) structure.
 - `type_3x3x3_confX_micmec.pickle` :
Basic output information, including the elastic stifness tensorn of the 3 x 3 x 3, micromechanical system for each of the eight configurations `X` visualized in Figure 7, as well as for the pure **fcu** (`X=0`) and the pure **reo** (`X=9`) structure.
- `output_configurations.py` :
Python script to extract from each of the aforementioned files the bulk modulus, Young moduli, and elastic stiffness tensors.
-  `output_configurations.json` :
Output file generated by the previous Python script.

### UiO-66/scripts_figures

This subfolder contains the scripts to generate Figures 8 and 9 from the above json output file:
- `fig_bulk_modulus` :
Python script to extract and visualize the bulk modulus of Figure 8a.
- `fig_young_modulus` :
Python script to extract and visualize the range of Young moduli of Figure 8b.
- `fig_dir_young_modulus` :
Python script to extract and visualize the Young modulus in 3D, as in Figure 9.
