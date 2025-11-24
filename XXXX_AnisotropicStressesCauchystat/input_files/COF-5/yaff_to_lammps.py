# -*- coding: utf-8 -*-
# YAFF is yet another force-field code.
# Copyright (C) 2011 Toon Verstraelen <Toon.Verstraelen@UGent.be>,
# Louis Vanduyfhuys <Louis.Vanduyfhuys@UGent.be>, Center for Molecular Modeling
# (CMM), Ghent University, Ghent, Belgium; all rights reserved unless otherwise
# stated.
#
# This file is part of YAFF.
#
# YAFF is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# YAFF is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#
# ----------------------------------------------------------------------

import os
#YAFF imports:
from yaff.external.lammpsio import write_lammps_table, ff2lammps
from yaff import System, ForceField
from molmod.units import angstrom

def main():
    input_dir = './yaff_input/'
    nx,ny,nz = 4,4,4
    # Generate supercell system
    system = System.from_file(input_dir + 'init.chk').supercell(nx,ny,nz)
    # Tabulate vdW interactions
    if not os.path.isfile('lammps.table'):
        ff = ForceField.generate(system, [input_dir + 'pars.txt'], rcut=15.0*angstrom)
        write_lammps_table(ff, fn='lammps.table', rmin=0.50*angstrom, nrows=2500, unit_style='real') #input: yaff forcefield, file name to write too, 
    # Write the LAMMPS input files
    ff2lammps(system, input_dir + 'pars.txt', ".", rcut=15.0*angstrom, tailcorrections=True,
        tabulated=True, unit_style='real')
    # Adapt the sampling options, which are defined in the last 5 lines
    # of lammps.in
    with open('lammps.in','r') as f:
        lines = f.readlines()
    with open('lammps.in','w') as f:
        # Write all lines except the last 10
        for line in lines[:-16]:
            # Insert 'pair_modify tail yes' before Output settings
            f.write(line)
            if 'pair_style' in line:
                f.write('pair_modify tail yes\n')
if __name__=='__main__':
    main()