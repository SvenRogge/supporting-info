


# Table of contents

This GitHub repository contains the input and raw output data accompanying the manuscript

**Dynamic layer stacking in 2D covalent organic frameworks and its suppression towards increased crystallinity**

by Sander Borgmans, XX, YY, and ZZ.


The data presented here is licensed under the CC BY-SA 4.0 international license, a copy of which can be found [here](https://creativecommons.org/licenses/by-sa/4.0/). Under this license, you can copy and redistribute the material in any medium or format as long as you give appropriate credit, provide a link to the license, and indicate if changes were made.

Additional information concerning the data is available upon request from the authors. Please send a mail to Sven.Rogge@UGent.be for more information.

## Software
*OGRe*
The free energy surfaces were obtained using OGRe, a Python package for optimal umbrella sampling grid refinement in an arbitrary number of dimensions. More information about OGRe can be found in [this manuscript](https://pubs.acs.org/doi/10.1021/acs.jctc.3c01028) and [this GitHub repository](https://github.com/SanderBorgmans/OGRe).

*Yaff and LAMMPS*
All molecular dynamics simulations were performed with Yaff, which was interfaced with LAMMPS to efficiently calculate the long-range interactions.

*Pre- and post-processing scripts*
The pre- and post-processing Python scripts are written in Python 3.12.


## Archive

This archive consists of five folders: one for each of the five COFs discussed in the manuscript.

Each of these folders contain the following files:

* `custom_cv`: a custom Python script that defines the two collective variables (CVs) defining the layer offset, as defined in the manuscript.
* `fes.dat`: the numerical free energy surface, containing three columns. These three columns correspond to the first CV (in Å), the second CV (in Å), and the resulting free enthalpy (in kJ mol<sup>-1</sup>) at that point. The files contain the numerical data that, once restricted to the Wigner-Seitz cell, correspond to the free enthalpy surfaces provided in Figures 2, 3, and 4 in the manuscript.
* `init.chk`: a Yaff input file containing the initial COF structure used in the umbrella sampling simulations
* `layer00.txt`, `layer01.txt`, and, for TTI-COF, `layer02.txt`: these so-called Layer files are generated sequentially by the OGRe protocol and contain information about the grid used to perform the umbrella sampling simulations. Each subsequent layer iteratively improves on the previous one. For TTI-COF, three layers were used, while two layers sufficed for all other COFs.
* `pars_XX.txt`: the force field parameters for the given structure. For COF-5, these are divided into the covalent (`pars_cov.txt`), the electrostatic (`pars_ei.txt`), and the van der Waals (`pars_mm3.txt`) parameters. For all other COFs, all parameters are provided in the file `pars_hcb_*.txt`, where the wildcard corresponds to the name of the COF structure according to [the ReDD-COFFEE nomenclature](https://pubs.rsc.org/en/content/articlelanding/2023/ta/d3ta00470h).