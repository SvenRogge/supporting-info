# Table of contents

This GitHub repository contains the input and raw output data accompanying the manuscript

**MOFs for long-term gas storage: exploiting kinetic trapping in ZIF-8 for on-demand and stimuli-controlled gas release**

by Karsten Heinz, Sven M.J. Rogge, Andreas Kalytta-Mewes, Dirk Volkmer, and Hana Bunzen.

This work was published in *Inorg. Chem. Front.*, **2023**, *10* (16): 4763-4772 DOI: [10.1039/D3QI01007D](https://dx.doi.org/10.1039/D3QI01007D)).

The data presented here is licensed under the CC BY-SA 4.0 international license, a copy of which can be found [here](https://creativecommons.org/licenses/by-sa/4.0/). Under this license, you can copy and redistribute the material in any medium or format as long as you give appropriate credit, provide a link to the license, and indicate if changes were made.

Additional information concerning this computational data is available upon request from the authors. Please send a mail to Sven.Rogge@UGent.be for more information.


## Software
*Gaussian*
The SF<sub>6</sub> molecule has been relaxed using the Becker 3-parameters hybrid exchange functional in combination with the Lee-Yang-Parr correlation (B3LYP) functional. The split-valence double-zeta 6-311G(d,p) basis set has been employed to describe the orbitals of the different atoms. All these methods have been adopted as implemented in the Gaussian 16 software suite, using the module g16_A.03-intel-2019a. The very tight convergence criteria were adopted for the optimization of the geometries.

*QuickFF and Horton*
To generate a force field for the SF<sub>6</sub> molecule, the HORTON module 2.1.1-intel-2020a-Python-2.7.18 and the QuickFF v2 module 2.2.7-intel-2020a-Python-3.8.2 were used.

*Yaff and LAMMPS*
All molecular dynamics simulations were performed with Yaff, which was interfaced with LAMMPS to efficiently calculate the long-range interactions. The  7Aug2019-foss-2019b-Python-3.7.4-kokkos LAMMPS module was used to this end. The probability histograms of the US simulations were converted to a free energy profile using a desktop-version of WHAM version 2.0.10.

*Pre- and post-processing scripts*
The pre- and post-processing Python scripts are written in Python 3.8.


## construct_ff_quickff_horton

### results
The folder `construct_ff_quickff_horton/results` contains the SF<sub>6</sub> force field (divided into a covalent, and electrostatic, and a van der Waals part) as well as the QuickFF System file for this molecule.

### input
The folder `construct_ff_quickff_horton/input` contains the structure of the SF<sub>6</sub> molecule from which the optimization was started, and the QuickFF input script.

### opt
The folder `construct_ff_quickff_horton/opt` contains the job script to perform the optimization, as well as a formatted checkpoint file of the optimized structure.

### freq
The folder `construct_ff_quickff_horton/freq` contains the optimized structure for which the Hessian calculation was started, the job script to calculate this Hessian, and the formatted checkpoint file of the resulting Hessian.


### mbis
The folder `construct_ff_quickff_horton/mbis` contains the Gaussian formatted checkpoint file used to determine the atomic charges, the job script to perform the MBIS calculation, and the output of this calculation.


### qff
The folder `construct_ff_quickff_horton/qff` contains the job script used to derive the QuickFF force field as well as the System chk files containing the different atom types.

## NPT_yaff_lammps

### results
This folder contains, for each loading between 1 and 6 SF<sub>6</sub> molecules in cage 1, a subfolder with initial System files. The `X` in `init_X.chk` indicates the value of the CV, in 0.1 Å. The prefix `m` is used to denote negative CV values.

### scripts
The `job.sh` file is an example of the files used to submit the molecular dynamics simulations to the HPC.

The `ymd.py` file  is an example of the input scripts used to start the molecular dynamics simulations.

The `test.py` file is used to generate a LAMMPS table of noncovalent interactions as well as to check that these interactions are consistent between Yaff and LAMMPS.

The `add_SF6` file is used to add a predefined number of SF<sub>6</sub> molecules to an empty ZIF-8 System object so that they don't overlap.

The `extract_first_snapshot.py` script is used to extract, from the here performed NPT runs, the first snapshot for which the observed CV value is close enough to a predefined CV value.

The `move_along_CV.py` script is used in case the above script does not yield an appropriate initial structure for a given CV value. In that case, this script extracts from the previous NPT simulation that structure for which the observed CV value is as close as possible to the required CV value, and further moves the SF<sub>6</sub> molecule as a rigid body along the normal that defines the CV.


### input
The `pars.txt` and the `ZIF8_4SF6.chk`/`ZIF8_6SF6.chk` files are the QuickFF force field and System objects used to perform a simulation with 6 SF<sub>6</sub> molecules inside the first cage of ZIF-8. These are constructed from the separate force fields (`pars_SF6.txt` and `pars_ZIF8.txt`) and from the separate System objects (`SF6.chk`and `ZIF8_empty.chk`), respectively.

The `lammps_smoothei2.table` file contains the LAMMPS tabulated data for the noncovalent interactions.


## US_NVT_yaff_LAMMPS

### results
This folder contains, for each of the six different loadings, the file `colvar_X_Y`, with the CV value as a function of the simulation time. In this nomenclature, `X` denotes the CV value (again in 0.1 Å) and `Y` is the bias strength, in kJ mol<sup>-1</sup>. In case `Y` is missing, the bias strength amounted to 25 kJ mol<sup>-1</sup>.

In addition, the file `CV_ZSF6_Y.csv` summarizes the CV values encountered when `Z` SF<sub>6</sub> molecules are present in cage 1, with `Y` again the bias strength as explained above. These data are also visualized.

Finally, `metadata.dat`, `free_energy.dat` and `free_energy.svg` are the metadata and output files from the WHAM analysis.

### scripts

The `job.sh` and `job_restart.sh` files are examples of the files used to submit the US molecular dynamics simulations to the HPC.

The files `ymd_us.py` and `ymd_us_restart.py` are needed to set up the US molecular dynamics simulations.

The file `calc_av.py` and `us_testoverlap.py` were used to verify whether sufficient overlap was present between adjacent umbrellas, both textually and visually.

The file `run_wham.py` was used to run the WHAM analysis.

The file `plot_all.py` was finally used to visualize all free energy profiles.
