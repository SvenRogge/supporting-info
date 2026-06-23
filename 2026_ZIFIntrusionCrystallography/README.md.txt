
# Table of contents

This GitHub repository contains the input data accompanying the manuscript

**Water clustering in hydrophobic zeolitic imidazolate frameworks under pressure**

by Titas Pramani, Jelto Neirynck, Daniel S. Parsons, Patrick W. Doheny, Hebin Jiang, Haixia Yin, Xiaochuan Liu, Ali Ashraf Siwji, Dominik Daisenberger, Mark R. Warren, David R. Allan, Sarah A. Barnett, Sven M. J. Rogge, Hamish H.-M. Yeung, and Yueting Sun.

It contains all the different input scripts for the *ab initio* optimisations and the machine-learned interatomic potential molecular dynamics used in this work.

The data presented here is licensed under the CC BY-SA 4.0 international license, a copy of which can be found [here](https://creativecommons.org/licenses/by-sa/4.0/). Under this license, you can copy and redistribute the material in any medium or format as long as you give appropriate credit, provide a link to the license, and indicate if changes were made.

Additional information concerning the data is available upon request from the authors. Please send a mail to Sven.Rogge@UGent.be for more information.

## Software
For all the DFT simulations in this work, _CP2K_ was used. _CP2K_ is available to download for free at [https://www.cp2k.org/download](https://www.cp2k.org/download). For the molecular dynamics simulations, the Langevin engine in _ASE_ was used in combination with the MACECalculator. ASE is available to download for free at [https://ase-lib.org/install.html](https://ase-lib.org/install.html). The ‘mace-mp-0b3-medium’ foundational model used in this work is available for free at https://github.com/ACEsuit/mace-foundations.


## Archive
This archive contains two folders: the folder *DFT_optimisations_CP2K* contains the required scripts and files for the two kinds of *ab initio* calculations: the geometry optimisations and the BSSE-corrected single point jobs, whereas the folder *MLIP_MD_ASE* contains the python script to generate the trajectories, and the starting structures used to generate the partial occupancy data.

### DFT_optimisations_CP2K
This folder contains four subfolders. 

The folder _CP2K_input_scripts_ contains the actual input scripts used in this manuscript for the three levels of constrained optimisation jobs. The file `cp2k_input_constrained_ZIF_O.inp` was used to optimise only the protons of the water molecules (stage 1 optimisations). In this script, the index “278” is included in &FIXED_ATOMS to fix the oxygen position of the water molecules during optimisations, as a water takes the “_H O H_” convention in this work. The file `cp2k_input_constrained_ZIF.inp` was used for the free water optimisations (stage 2 optimisations), and `cp2k_input_constrained_free.inp` was used for the fully flexible geometry optimisations (stage 3 optimisations). In addition, `cp2k_input_BSSE.inp` contains the input for the BSSE interaction energy calculations to compute the interaction energies. Currently, the second fragment in is left open, and needs to be filled in with the desired water cluster to compute the interaction energy. 

The folder _CP2K_resources_ contains the basis set (`basis_set`), dispersion correction (`dftd3`), and pseudopotentials (`potential`) required to run the *CP2K* jobs. 

The folder *input_structures* contains all the initial configurations from which the optimisations started, with for the single-molecule optimisations the notation `ZIF8_xmol_Oy_vz_θ_φ.xyz`, in which
* `x` represents the number of water molecules, in between 1 and 4, with 1 reported in the main text and 2-4 in the Supplementary Information;
* `Oy` represents one of the five adsorption sites O1 to O5;
* `z` varies between 0 and 2 and distinguishes between the three different crystallographic positions for the O1 and O4 sites;
* `θ` and `φ`: describes the rotation of the water molecules with respect to the Cartesian axis system (*x*, *y*, *z*) defined by the simulation box where the cell vector ***a*** is aligned along *x*, ***b*** along *y*, and ***c*** along *z*. Specifically, `θ` is the angle between the O-H bond of the first proton with respect to the `-z` axis, whereas `φ-21.6°` is the angle with respect to the `x` axis.

A similar but repeated notation was used for the optimisations of multiple water molecules.

Finally, the folder *minimum_energy_structures* contains the minimum-energy structure for each of the five adsorption sites when a single water molecule is present, and this for all three types of optimisations, as well as the minimum-energy structure from the stage 2 optimisations for two, three, and four water molecules.


### MLIP_MD_ASE
This folder contains the script `run_mace.py` and the subfolder `input_structures`.

The script `run_mace.py` was used to run the MLIP MD simulations using the "mace-mp-0b3-medium" model, available at https://github.com/ACEsuit/mace-foundations. The script will generate a single trajectory file, a checkpoint file, and a log file containing information on the temperature and potential energy. This script can be run by calling `python stability_run_mace.py [STRUCTURE_PATH] [TEMPERATURE] [OUTDIR]`, with `[STRUCTURE_PATH]` the path to the starting structure, `[TEMPERATURE]` the desired NVT temperature, and `[OUTDIR]` the desired folder for the trajectory.

The subfolder `input_structures` contains the initial structures for the MD simulations, all named `ZIF_8_x_y_vz`, where
* `x` represents the number of water molecules in the first cage;
* `y` represents the number of water molecules in the second cage;
* `z` varies between 0 and 2 and distinguishes between the three independent initial structures for the MD simulations per thermodynamic condition.