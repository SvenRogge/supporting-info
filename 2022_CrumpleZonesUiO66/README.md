
# Table of contents

This GitHub repository contains the input and raw output data accompanying the manuscript

**Absorbing stress via molecular crumple zones: Strain engineering flexibility into the rigid UiO-66 material**

by Sven M.J. Rogge, Sander Borgmans, and Veronique Van Speybroeck.

The data presented here is licensed under the CC BY-SA 4.0 international license, a copy of which can be found [here](https://creativecommons.org/licenses/by-sa/4.0/). Under this license, you can copy and redistribute the material in any medium or format as long as you give appropriate credit, provide a link to the license, and indicate if changes were made.

Additional information concerning the data is available upon request from the authors. Please send a mail to Sven.Rogge@UGent.be for more information.

## Software
*QuickFF*
To generate the periodic force fields for the six UiO-66 materials, QuickFF v2 module 2.2.7-intel-2020a-Python-3.8.2 was used.

*Yaff and LAMMPS*
All molecular dynamics simulations were performed with Yaff, which was interfaced with LAMMPS to efficiently calculate the long-range interactions. The  7Aug2019-foss-2019b-Python-3.7.4-kokkos LAMMPS module was used to this end.

*Pre- and post-processing scripts*
The pre- and post-processing Python scripts are written in Python 3.8.


## Archive

This archive consists of seven folders: one containing the relevant pre- and postprocessing scripts, and six containing the specific input and output files for each of the six UiO-66 materials. The latter six folders are named *CN**, where CN stands for coordination number. The wildcard corresponds to:
* *CN12*: The defect-free material of Fig. 1a in the manuscript, whose analysis is shown in Fig. 2 of the manuscript.
* *CN11-5*: The strain-engineered material of Fig. 1b in the manuscript containing a single linker vacancy, whose analysis is shown in Fig. S1 of the Supplemental Information.
* *CN11*: The strain-engineered material of Fig. 1c in the manuscript containing two linker vacancies, whose analysis is shown in Fig. 3 of the manuscript.
* *CN10-5a*: The strain-engineered material of the left panel of Fig. 1d in the manuscript containing three linker vacancies, whose analysis is shown in Fig. S2 of the Supplemental Information.
* *CN10-5b*: The strain-engineered material of the right panel of Fig. 1d in the manuscript containing three linker vacancies, whose analysis is shown in Fig. S3 of the Supplemental Information.
* *CN10*: The strain-engineered material of Fig. 1e in the manuscript containing four linker vacancies, whose analysis is shown in Fig. 4 of the manuscript.

The structure of each of these folders is discussed below.

### scripts

The subfolder **NPsT** contains the scripts used for the $(N, P, \bm \sigma_a = \bm 0, T)$ simulations. The `ymd.py` script performs the molecular dynamics simulation at a predefined pressure and temperature, and is started using `job.sh`. Afterwards, `check_volumes.py` is run to check whether all necessary volumes in a predefined range are encountered during the simulation. If so, `makefile.py` extracts snapshot at the specified volumes to create the pressure-versus-volume equations of state from the generated output files. Finally, `extrapolate_volume.py` is used as an alternative approach to generate initial structure files by extrapolating from a predefined volume, which is  used to test the robustness of the method in Section S2.

The subfolder **NVsT** contains the scripts used for the $(N, V, \bm \sigma_a = \bm 0, T)$ simulations. The `ymd.py` script performs the molecular dynamics simulation at a predefined temperature and is started using `job.sh`.  To create the strain fields visualized in panel a of each figure in the manuscript, `determine_Zr_clusters.py` is first run to determine the different groups of six zirconium atoms that form an inorganic building block, accounting for periodic boundary conditions. These clusters are then stored in `clusters.csv` (*vide infra*). Based on this and the MD output file, `locate_subcells.py` defines the different parallelepipeds $\{\mathbf{h}_{\mu \nu \kappa} (t) \}$, from which the strain fields are calculated and visualized using `calculate_strainfield.py`. It is important to first run this last script on the equilibrium structure, since this is necessary to define the reference structure $\mathbf{h}_\text{ref}$. The scripts `PvsV.py` and `PvsV.sh` are used to construct the pressure-versus-volume equations of state, which are then visualized using `plot_free_energy.py`. Finally, `extract_porosity.py` is used to extract the gravimetric accessible pore volume, which is visualized using `plot_porosity.py`.

#### CNXX

The subfolder **CNXX/ff** contains the periodic force field files used for the specific structure. The files `pars_cov.txt`, `pars_ei.txt`, and `pars_mm3.txt` contain the covalent, electrostatic, and van der Waals interactions, respectively.

The subfolder **CNXX/NPsT/input** contains the `init.chk` file, which describes the force-field optimized structure of the specified material. This is used as input structure for the $(N, P, \bm \sigma_a = \bm 0, T)$ simulations.

The subfolder **CNXX/NPsT/output** (which would largely coincide with the **CNXX/NVsT/input** subfolder, which is therefore omitted) contains the files `init_XXXX.chk`, which are the structure files extracted from the $(N, P, \bm \sigma_a = \bm 0, T)$ simulations at the volume `XXXX` Å$^3$ and which are used as input for the subsequent $(N, V, \bm \sigma_a = \bm 0, T)$ simulations. In case the robustness of the pressure-versus-volume equations of state are tested (see Section S2), this subfolder also contains files of the type `init_XXXX_extrapolated_from_YYYY.chk`. In this case, the structure at the volume `XXXX` Å$^3$ is not extracted directly from the $(N, P, \bm \sigma_a = \bm 0, T)$ simulation, but is rather extrapolated from the earlier extracted structure at the volume `YYYY` Å$^3$. Finally, for the CN10 case, this subfolder also contains `traj_1000MPa_strain.h5`, which contains the strain field of this structure at a pressure of 1000 MPa, discussed in Fig. 5 of the manuscript and provided as a Supplemental Video, in the HDF5 format.

The subfolder **CNXX/NVsT/output** contains `clusters.csv` as output from the `determine_Zr_clusters.py` script discussed above. Each different row of six atomic indices corresponds to a different zirconium cluster. The file `porosity.csv` is the output from the script `extract_porosity.py` and contains, from left to right, the supercell volume (in Å$^3$, to be divided by eight to get the conventional unit cell volume), the density of the material, and the accessible pore volume (in Å$^3$, volume fraction, and cm$^3$g$^{-1}$). Only the latter is reported in the manuscript, as all other accessible pore volume metrics sketch a similar picture. The file `PvsV_eq100ps_av900ps.csv` contains the supercell volume (first column, in Å$^3$, to be divided by eight to get the conventional unit cell volume) and the average pressure at that volume (second column, in MPa), as extracted via `PvsV.py`. In case different initial structures are used (via interpolation), additional  files of the form `PvsV_eq100ps_av900ps_extrapolated_from_YYYY.csv` are present, where `YYYY` adopts the same meaning as above. Finally, `traj_XXXX_strain.h5` contains the strain field at the volume `XXXX` Å$^3$, discussed in panel a of Figs. 2, 3, 4, S1, S2, and S3, and provided as Supplemental Videos.
