# Table of contents

This GitHub repository contains the input and raw output data accompanying the manuscript

**Unraveling the Thermodynamic Criteria for Size-Dependent Spontaneous Phase Separation in Soft Porous Crystals**

by Sven M.J. Rogge, Michel Waroquier, and Veronique Van Speybroeck.

This work was published in *Nat. Commun.*, **2019**, *10*: 4842 (DOI: [10.1038/s41467-019-12754-w](dx.doi.org/10.1038/s41467-019-12754-w)).

The data presented here is licensed under the CC BY-SA 4.0 international license, a copy of which can be found [here](https://creativecommons.org/licenses/by-sa/4.0/). Under this license, you can copy and redistribute the material in any medium or format as long as you give appropriate credit, provide a link to the license, and indicate if changes were made.

Additional information concerning the data is available upon request from the authors. Please send a mail to Sven.Rogge@UGent.be for more information.


## MIL-53(Al)

The MIL-53(Al) folder contains three subfolders:

1. *ff*: contains the force field and the initial geometry used to perform the simulations
2. *NPsT*: contains the python scripts for the $(N,P,  \bm \sigma_a =\bm 0, T)$  simulations
3. *NVsT*: contains the python scripts for the $(N,V,  \bm \sigma_a =\bm 0, T)$  simulations, as well as the intermediate and final results of the pressure equations of state for each of the supercells. These files contain two columns: the first column is the supercell volume (in ångström), the second column the averaged internal pressure (in MPa)


## DMOF-1(Zn)

The DMOF-1(Zn) folder contains three subfolders:

1. *ff*: contains the force field and the initial geometry used to perform the simulations
2. *NPsT*: contains the python scripts for the $(N,P,  \bm \sigma_a =\bm 0, T)$  simulations
3. *NVsT*: contains the python scripts for the $(N,V,  \bm \sigma_a =\bm 0, T)$  simulations, as well as the intermediate and final results of the pressure equations of state. These files contain two columns: the first column is the supercell volume (in ångström), the second column the averaged internal pressure (in MPa)


## MIL-53(Al)-F

The MIL-53(Al)-F folder contains three subfolders:

1. *ff*: contains the force field and the initial geometry used to perform the simulations
2. *NPsT*: contains the python scripts for the $(N,P,  \bm \sigma_a =\bm 0, T)$  simulations
3. *NVsT*: contains the python scripts for the $(N,V,  \bm \sigma_a =\bm 0, T)$  simulations, as well as the intermediate and final results of the pressure equations of state at 100 K, 300 K, and 500 K. These files contain two columns: the first column is the supercell volume (in ångström), the second column the averaged internal pressure (in MPa)

## CoBDP

The CoBDP folder contains three subfolders:

1. *ff*: contains the force field and the initial geometry used to perform the simulations
2. *NPsT*: contains the python scripts for the $(N,P,  \bm \sigma_a =\bm 0, T)$  simulations
3. *NVsT*: contains the python scripts for the $(N,V,  \bm \sigma_a =\bm 0, T)$  simulations, as well as the intermediate and final results of the pressure equations of state at methane loadings of 0, 2, and 4 molecules per conventional unit cell. These files contain two columns: the first column is the supercell volume (in ångström), the second column the averaged internal pressure (in MPa)

## CIFs

The CIFs folder contains the crystallographic information files for the different phase coexistence regions and pure phases considered in the manuscript (see Supplementary Note 9). The folder contains six subfolders:

1. *MIL53Al_424*: contains the CIFs for the 4x2x4 cell of MIL-53(Al) at 300 K
2. *MIL53Al_626*: contains the CIFs for the 6x2x6 cell of MIL-53(Al) at 300 K
3. *MIL53Al_828*: contains the CIFs for the 8x2x8 cell of MIL-53(Al) at 300 K
4. *DMOF1Zn*: contains the CIFs for the 8x2x8 cell of DMOF-1(Zn) at 300 K
5. *MIL53AlF_100K*: contains the CIFs for the 8x2x8 cell of MIL-53(Al)-F at 100 K
6. *CoBDP_0CH4*: contains the CIFs for the empty 8x2x8 cell of CoBDP at 300 K


