# Table of contents

This GitHub repository contains the input and raw output data accompanying the manuscript

**"High rate nanofluidic energy absorption in porous zeolitic frameworks"**
by Yueting Sun, Sven M. J. Rogge, Aran Lamaire, Steven Vandenbrande, Jelle Wieme, Clive R. Siviour, Veronique Van Speybroeck, and Jin-Chong Tan.

The data presented here is licensed under the CC BY-SA 4.0 international license, a copy of which can be found [here](https://creativecommons.org/licenses/by-sa/4.0/). Under this license, you can copy and redistribute the material in any medium or format as long as you give appropriate credit, provide a link to the license, and indicate if changes were made.

Additional information concerning the data is available upon request from the authors. Please send a mail to Sven.Rogge@UGent.be for more information.


## Software
*VASP*
The periodic ZIF-8 material was relaxed with VASP (module 5.4.1-intel-2015b-mt-vaspsol-20150914) using the projector-augmented wave (PAW) method. The PBE exchange-correlation functional was combined with the DFT-D3 dispersion scheme using Becke-Johnson damping. The recommended GW PBE PAW potentials were employed for all elements and functionals (v5.4). For the zinc atoms, the 3s, 3p, 3d, and 4s electrons were included explicitly. For the carbon and nitrogen atoms, the 2s and 2p electrons were considered as valence electrons. For the hydrogen atoms, the 1s electron was treated as a valence electron. These calculations were performed with a plane-wave kinetic-energy cut-off of 800 eV and using Gaussian smearing with a smearing width of 0.05 eV. Projection operators were evaluated in reciprocal space. A Γ-point k-grid was used for all volumes. The real-space FFT grid was used to describe wave vectors up to two times the maximum wave vector present in the basis set. An augmentation grid that is twice as large was used to avoid wrap-around errors in order to obtain accurate forces. The electronic (ionic) convergence criterion was set to 10<sup>-9</sup> (10<sup>-8</sup>) eV. Subsequently, the dynamical matrix was determined using 0.015 &#8491; displacements for all atomic coordinates with respect to the equilibrium structure.

*GPAW and HORTON*
The GPAW module 0.11.0.13004-intel-2015b-Python-2.7.10 and the corresponding GPAW-setups module 0.9.11271-linux-x86_64 have been employed to determine the charges using the PBE exchange-correlation functional and a grid spacing of 0.15 &#8491;. The resulting charges were then extracted using the MBIS scheme as implemented in HORTON. 

*QuickFF*
From the dynamical matrix and the charges determined above, a flexible and _ab initio_-based force field for the empty ZIF-8 structure was derived. The covalent part of the force field, which contains diagonal terms that describe bonds, bends, out-of-plane distances, and torsion angles, as well as cross terms was derived using the in-house QuickFF software package (version 2.2.0). The van der Waals interactions in these force field models were added *a posteriori* based on the MM3 parameters.


*CP2K*
In order to probe the influence of temperature on the ZIF-8 swing angle, an additional set of _ab initio_ MD simulations were performed at temperatures of 100 K, 200 K, and 300 K and at a pressure of 0 MPa using the CP2K software package (module 5.1-intel-2018a). In these calculations, the PBE-D3(BJ) level of theory was used in combination with Gaussian TZVP-MOLOPT basis sets, a plane wave basis set with a cut-off of 800 Ry and a relative cut-off of 60 Ry, and Goedecker-Teter-Hutter (GTH) pseudopotentials. The temperature of the simulations was controlled with a Nosé-Hoover chain thermostat consisting of three beads and with a time constant of 0.1 ps. The pressure was controlled with a Martyna-Tobias-Klein barostat  with a time constant of 1 ps. The MD time step was set to 0.5 fs. The total simulation time for the AIMD simulations comprised 11 ps, of which the first ps was discarded for equilibration.

*RASPA2*
To determine the water saturation and obtain representative snapshots of water-filled ZIF-8 used in the force-field MD simulations, grand canonical Monte Carlo (GCMC) simulations in the NVT ensemble were performed using the RASPA2 software package (module 2.0.3-vtkrestart-intel-2018a). All simulations were performed at 298 K and various pressures. During these MC simulations, the framework and the internal coordinates of the water molecules were kept rigid. The Lennard-Jones interactions were truncated at 12 Å with analytical tail corrections to correct for the finite cut-off. Electrostatic interactions were treated using the Ewald summation method.

*Yaff*
All force-field-based molecular dynamics simulations were performed with Yaff 1.4.2, using the module 1.4.2-intel-2017b-Python-3.6.3. The MD simulations were run for a total simulation time of 3 ns (inhomogeneous water distribution) or 5 ns (homogeneous water distribution). During these simulations, the temperature was controlled to be on average 300 K using a Nosé-Hoover chain thermostat consisting of three beads and with a time constant of 0.1 ps. The pressure was controlled with a Martyna-Tobias-Klein barostat with a time constant of 1 ps. The integration time step was limited to 0.5 fs to ensure energy conservation when using the velocity Verlet scheme. The long-range van der Waals interactions were cut off at a radius of 12 Å, which was compensated by tail corrections. The electrostatic interactions were efficiently calculated using an Ewald summation with a real-space cutoff of 12 Å, a splitting parameter α of 0.213 Å<sup>-1</sup>, and a reciprocal space cutoff of 0.32 Å<sup>-1</sup>.

*Pre- and post-processing scripts*
The pre- and post-processing Python scripts are written in Python 2.7 or 3.7. To convert the trajectory files to crystallographic information files (CIFs), CON3F was used (module version_CIF2-intel-2015a-Python-2.7.10).


## construct_ff_quickff_vasp_horton

### input
This folder contains the VASP input files (`CONTCAR`, `INCAR`, and `KPOINTS`), the GPAW starting structure (`gpaw_cluster.xyz`), and the QuickFF configuration file  (`quickff_input.txt`).

### results
This folder contains the final results of the force field derivation. The `OUTCAR` file contains the VASP output file,  while the `pars.txt` and the `system_opt.chk` files are the force field an QuickFF system file after force-field optimization. 

### scripts
This folder contains the job scripts to start the VASP, GPAW, HORTON, and QuickFF simulations (`run_vasp.sh`, `run_gpaw.sh`, `run_mbis_denspart.sh`, and `run_quickff.sh`, respectively). In addition, `run_yaff_opt.py` is the Python script used to carry out the initial Yaff optimization of the ZIF-8 structure.

## gcmc_raspa2

### input
This folder contains four subfolders.

The **conversions** subfolder contains text files summarizing the most important conversion rules used to report the uptake of guest molecules for each type of guest molecule and both for the AP and HP phases of ZIF-8.

The **ff** subfolder contains the ZIF-8 force field file `pars.txt` derived above, as well as the definition of the force field and force field mixing rules used in RASPA2 (`force_field.def` and `force_field_mixing_rules.def`).

The **systems** subfolder contains for each of the guest molecules the RASPA2 input files (`*.def`) as well as the Yaff system files (`*.chk`), the latter used to define the force field parameters in the MD simulations below. In addition, it contains the ZIF-8 CIFs used as input for the GCMC simulations (`ZIF8_AP.cif` and `ZIF8_HP.cif`, corresponding to the CCDC files with codes TUDHUW and TUDJOS) and the accompanying Yaff system files (`*.chk`). Finally, `pseudo_atoms.def` contains all atoms that are not explicitly moved during the GCMC simulations (the framework and ghost atoms).

The **template_scripts** subfolder contains template input scripts to construct the RASPA2 simulation grid (`makegrid.input`), to perform the insertion simulations (`insertion.input`), and to start and restart the GCMC simulations (`simulation.input` and `restart.input`). Furthermore, `temperatures.dat` and `pressures.dat` contain a list of temperatures and pressures for which simulations have been carried out (not all temperatures and pressures are directly relevant for this manuscript).

### results

This folder contains three subfolders as well as the file `summary.txt` that summarizes all simulations performed for this work and the resulting uptake.

The **plot** subfolder contains three types of plots. The subfolders **ZIF8_AP_111**,  **ZIF8_AP_222**,  **ZIF8_HP_111**, and  **ZIF8_HP_222** visualizes for each of the relevant amount of guest molecules the convergence of the simulations. Note that for the supercells also the guest loading is multiplied by eight.  The **ZIF8_comparison_111** and **ZIF8_comparison_222** folders compare the converged densities between the AP and HP structures. For these six subfolders, data that is symmetrized over the different unit cells is also available in the subfolders with the suffix ***_symm**.

### scripts

This folder contains three subfolders.

The **ff** subfolder contains the implementation of the ghost atoms (`ghostatoms.py`) and tail corrections (`tailcorr.py`) used during the simulation. Furthermore, `make_guests.py` is used to generate a single guest molecule with the correct force field atom types resulting in the `*.chk` files in the input folder. Finally,  `write_raspa_input.py` is used to automatize the generation of the RASPA2 `*.def` files, calling `cif.py` and `medff.py`.

The **start_simulations** subfolder contain a job script for each of the three job types (initalization of the simulation grid (`run_makegrid.sh`), insertion simulations (`run_insertion.sh`) and GCMC simulations (`run_gcmc.sh`)) as well as three overarching job scripts to follow up these simulations (`job_*.sh`). The file `insertion.py` is the custom-made Python script to insert a given amount of guest molecules in the system.

The **post_analysis** subfolders contains all post-processing files. The file `raspa_grids.py` is used to convert the RASPA2 output files to `*.cube` files, which are processed by `plot_densities.py` to generate density profiles. The files `write_results.py` and `write_results.sh` are used to extract and summarize relevant uptake data from the RASPA2 output files, which are then used by `write_isotherms.sh` to write a summarizing text file that can be used to visualize the corresponding isotherm. The `process_loading.py` script outputs the guest loading during a RASPA2 simulation. Finally, `write_snapshots.py` is used to extract from a RASPA2 GCMC simulation a snapshot with a predetermined amount of guest molecules.

## NPT_abinitio_cp2k

### input
This folder contains an example CP2K input script (`input.inp`, used for the 100 K simulation), as well as the file containing the ZIF-8 structure used to initialize the CP2K simulation (`zif8_opt.xyz`).

### results
The `CP2K_angles_X00K.dat` files contain the time evolution of all dihedral swing angles at the simulation temperatures of 100 K, 200 K, and 300 K. Each column corresponds to the time evolution of one specific dihedral swing angle.

### scripts
This folder contains a typical CP2K job script (`job.sh`), as well as the file `determine_rings.py`, which processes the CP2K output data and outputs the time evolution of the dihedral swing angles during the simulation.


## NPT_forcefield_yaff

### input

This folder contains three subfolders as well as the files `pars_ZIF8.txt` and `pars_ZIF8_H2O.txt` that contain the QuickFF derived force fields with and without the terms pertaining to the TIP4P water. 

The **ZIF8_empty** subfolder contains the initial system file used to start the Yaff simulations of the empty framework (`init.chk`).

The **water_homogeneous** subfolder contains the initial system files with various amounts of water molecules inside ZIF-8 (between 0 and 80 molecules per 1x1x1 unit cell) , used to start the Yaff simulations for the homogeneous water distributions (`init_XXX.chk`).

The **water_gradient** subfolder contains the initial system file with the inhomogeneous water distribution (42 water molecules in a given pore of a 1x1x2 supercell, all other pores empty) considered in the manuscript (`init.chk`).

### results

This folder contains two subfolders. Note that the results for the inhomogeneous water distribution are not available here since their size is too large; but they are available upon request.

The **ZIF8_empty** subfolder contains the time evolution of the empty ZIF-8 structure under various pressures in  Å<sup>3</sup> (`volume_ZIF8_empty_XXXMPa.txt`).

The **water_homogeneous** subfolder contains, for each pressure and water loading at 300 K, the 20 CIFs obtained by splitting the total simulation of 5 ns in 20 subsimulations of 250 ps each and averaging the structure within each of these 250 ps subsimulations (`traj_XXX_XX.cif`). 

### scripts

This folder contains four subfolders.

The **ff** subfolder contains the implementation of the ghost atoms present in the TIP4P force field (`ghostatoms.py`) and the analytical tail corrections used during the simulations (`tailcorr.py`). Both files are called by the various `ymd.py` and `ymd_restart.py` scripts below.

The **ZIF8_empty** subfolder contains a representative Python script to run the Yaff simulation (`ymd.py`) as well as the accompanying job script (`job.sh`). Furthermore, it contains the script `read_volume.py`, which extracts the cell volume from the Yaff-generated `*.h5` files.

The **water_homogeneous** subfolder contains representative Python scripts to start and restart the Yaff simulations (`ymd.py` and `ymd_restart.py`) as well as the accompanying job scripts (`job.sh` and `job_restart.sh`). Furthermore, it contains scripts to generate (`generate_*.py`) and plot (`plot_*.py`) the dihedral swing angle distribution (`*_dihedrals.py`), the mean-squared displacement (`*_msd.py`), the radial distribution function with respect to the central pore center (`*_rdf.py`), the XRD patterns, generated from the CIFs (`generate_xrd.py`), and the number of hopping events (`plot_number_hoppings.py`).

The **water_gradient** subfolder contains representative Python scripts to start and restart the Yaff simulations (`ymd.py` and `ymd_restart.py`) as well as the accompanying job scripts (`job.sh` and `job_restart.sh`). Furthermore, it contains the script to generate the mean-squared displacement of water in the 1x1x2 ZIF-8 supercell (`generate_msd.py`).
