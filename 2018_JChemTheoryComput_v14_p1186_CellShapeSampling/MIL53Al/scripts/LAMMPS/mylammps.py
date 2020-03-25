# -*- coding: utf-8 -*-
# YAFF is yet another force-field code
# Copyright (C) 2011 - 2013 Toon Verstraelen <Toon.Verstraelen@UGent.be>,
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
#--
'''Lammps interface

   The mylammps module contains an interface to outsource the computation of
   non-covalent forces to the Lammps library.
'''

import numpy as np, os
import matplotlib.pyplot as pt
import ctypes
from scipy.special import erf
from mpi4py import MPI

from molmod.units import pascal, angstrom, kjmol, kcalmol
from molmod.constants import boltzmann

from yaff.log import timer
from yaff.pes.ff import ForcePart
from yaff.pes.ext import Cell, nlist_recompute, PairPotEI, Switch3
from yaff.pes import Scalings, ForcePartPair
from yaff.pes.nlist import NeighborList

from liblammps import lammps

#TODO 
#Included BondedNeighborList
__all__=[
    'ForcePartLammps', 'write_lammps_table', 'write_lammps_data', 'BondedNeighborList'
]


neigh_dtype = [
    ('a', int), ('b', int), ('d', float),        # a & b are atom indexes, d is the distance
    ('dx', float), ('dy', float), ('dz', float), # relative vector (includes cell vectors of image cell)
    ('r0', int), ('r1', int), ('r2', int)        # position of image cell.
]

class BondedNeighborList(NeighborList):
    '''A neighbor list that is intended for near-neighbor interactions. The
       pairs in the list are never updated, only distances are recomputed.
    '''
    def __init__(self, system, selected=None):
        '''
           **Arguments:**

           system
                A System instance.

           **Optional arguments:**

           selected
                A NumPy array [nneigh x 2] containing all pairs that should
                be considered.
                Default: All 1-2, 1-3 and 1-4 pairs included
        '''
        self.system = system
        if selected is None:
            bonded = []
            # Add all 1-2 neighbors
            for bond in system.bonds:
                bonded.append(bond)
            # Add all 1-3 neighbors
            for angle in system.iter_angles():
                bonded.append(np.array([angle[0],angle[2]]))
            # Add all 1-4 neighbors
            for dihed in system.iter_dihedrals():
                pass#bonded.append(np.array([dihed[0],dihed[3]]))
            for i in xrange(len(bonded)):
                bonded[i] = [np.amax(bonded[i]),np.amin(bonded[i])]
            bonded = np.asarray(bonded)
            # Only retain unique pairs
            selected = np.array([np.array(x) for x in set(tuple(x) for x in bonded)])
        self.nneigh = selected.shape[0]
        self.neighs = np.empty((self.nneigh),dtype=neigh_dtype)
        for ibond, bond in enumerate(selected):
            self.neighs[ibond]['a'] = np.amax(bond)
            self.neighs[ibond]['b'] = np.amin(bond)
            self.neighs[ibond]['r0'] = 0
            self.neighs[ibond]['r1'] = 0
            self.neighs[ibond]['r2'] = 0

    def request_rcut(self, rcut):
        # Nothing to do...
        pass

    def update_rmax(self):
        # Nothing to do...
        pass

    def update(self):
        # Simply recompute distances, no need to rebuild
        nlist_recompute(self.system.pos, self.system.pos, self.system.cell, self.neighs[:self.nneigh])
 
class ForcePartLammps(ForcePart):
    '''Energies obtained from Lammps.'''
    def __init__(self, system, pppm_accuracy=1e-5, fn_log='lammps.log', scalings=np.zeros(6),
                    fn_system='lammps.data', fn_table='lammps.table', triclinic=True, comm=None):
        '''
           **Arguments:**

           system
                An instance of the ``System`` class.

           **Optional Arguments:**

           scalings
                Numpy array [6x1] containing the scaling factors for 1-2, 1-3,
                1-4 Lennard-Jones and 1-2, 1-3, 1-4 electrostatic interactions.
                Default: [0.0,0.0,0.0,0.0,0.0,0.0]

           pppm_accuracy
                Desired relative error in electrostatic forces
                Default: 1e-5

           fn_log
                Filename where LAMMPS output is stored.
                Default: 'lammps.log'

           fn_system
                Filename of file containing system information, can be written
                using the ```write_lammps_data``` method.
                Default: lammps.data

           fn_table
                Filename of file containing tabulated non-bonded potential
                without charges, can be written using the
                ```write_lammps_table``` method.
                Default: lammps.table

           triclinic
                Boolean, specify whether a triclinic cell will be used during
                the simulation. If the cell is orthogonal, set it to False
                as LAMMPS should run slightly faster.
                Default: True

           comm
                MPI communicator, required if LAMMPS should run in parallel
        '''
        if system.cell.nvec != 3:
            raise ValueError('The system must be 3d periodic for Lammps calculations.')
        if not os.path.isfile(fn_system):
            raise ValueError('Could not read file %s' % fn_system)
        if not os.path.isfile(fn_table):
            raise ValueError('Could not read file %s' % fn_table)
        ForcePart.__init__(self, 'lammps', system)
        self.system = system
        self.comm = comm
        self.triclinic = triclinic
        self.setup_lammps(fn_system, fn_table, pppm_accuracy, fn_log, scalings)
        # LAMMPS needs cell vectors (ax,0,0), (bx,by,0) and (cx,cy,cz)
        # This means we need to perform a rotation to switch between Yaff and
        # LAMMPS coordinates. All information about this rotation is stored
        # in the variables defined below
        self.rvecs = np.eye(3)
        self.cell = Cell(self.rvecs)
        self.rot = np.zeros((3,3))

    def setup_lammps(self, fn_system, fn_table, pppm_accuracy, fn_log, scalings):
        '''
        Pass all commands that would normally appear in the LAMMPS input file
        to our instance of LAMMPS.
        '''
        self.lammps = lammps(name='mpi', comm=self.comm, cmdargs=["-screen","none","-log",fn_log])
#        self.lammps = lammps(name='mpi', comm=self.comm)
        nffa = self.system.ffatypes.shape[0]
        self.lammps.command("units electron")
        self.lammps.command("atom_style full")
        self.lammps.command("atom_modify map array")
        self.lammps.command("box tilt large")
        self.lammps.command("read_data %s"%fn_system)
        self.lammps.command("mass * 1.0")
        self.lammps.command("bond_style none")
        if self.system.charges is not None:
            self.lammps.command("pair_style hybrid/overlay coul/long 28.4 table spline 5000")
            self.lammps.command("pair_coeff * * coul/long")
            self.lammps.command("kspace_style pppm %f" % pppm_accuracy)
        else:
            self.lammps.command("pair_style table spline 2000")
#        # Electrostatics only        
#        self.lammps.command("pair_style coul/long 24.0")
#        self.lammps.command("pair_coeff * *")
#        self.lammps.command("kspace_style pppm %f" % pppm_accuracy)
        for i in xrange(nffa):
            ffai = self.system.ffatypes[i]
            for j in xrange(i,nffa):
                ffaj = self.system.ffatypes[j]
                if ffai>ffaj:
                    name = '%s-%s' % (ffai,ffaj)
                else:
                    name = '%s-%s' % (ffaj,ffai)
                self.lammps.command("pair_coeff %d %d table %s %s" % (i+1,j+1,fn_table,name))
#                self.lammps.command("pair_coeff %d %d %s %s" % (i+1,j+1,fn_table,name))
#                self.lammps.command("pair_coeff %d %d %s %03d-%03d" % (i+1,j+1,fn_table,i,j))
        if self.system.charges is not None:
            self.lammps.command("special_bonds lj %f %f %f coul %f %f %f" %
                 (scalings[0],scalings[1],scalings[2],scalings[3],scalings[4],scalings[5]))
        else:
            self.lammps.command("special_bonds lj %f %f %f" %
                 (scalings[0],scalings[1],scalings[2]))
        self.lammps.command("neighbor 0.0 bin")
        self.lammps.command("neigh_modify delay 0 every 1 check no")
        self.lammps.command("variable eng equal pe")
        self.lammps.command("compute virial all pressure NULL virial")
        self.lammps.command("fix 1 all nve")

    def update_rot(self):
        # Compute the transformation to go from Yaff to LAMMPS coordinates,
        # based on current Yaff and LAMMPS cell vectors
        self.rot[:] = 0.0
        A = self.system.cell.rvecs[0,:]
        B = self.system.cell.rvecs[1,:]
        C = self.system.cell.rvecs[2,:]
        self.rot[0,:] = np.cross(B,C)
        self.rot[1,:] = np.cross(C,A)
        self.rot[2,:] = np.cross(A,B)
        self.rot = np.dot(self.rvecs.transpose(),self.rot)/self.system.cell.volume

    def update_pos(self, pos):
        '''
        Update the LAMMPS positions based on the coordinates from Yaff
        '''
        # Perform the rotation
        pos[:] = np.einsum('ij,kj', pos, self.rot)
        # TODO: check if mic is necessary or not
#        for i in xrange(self.system.natom):
#            self.cell.mic(pos[i])
#        x = self.lammps.gather_atoms("x",1,3)
#        for i in xrange(3*self.system.natom):
#            x[i] = pos[i/3,i%3]
        #self.lammps.scatter_atoms("x",1,3,x)
        self.lammps.scatter_atoms("x",1,3,ctypes.c_void_p(pos.ctypes.data))

    def update_rvecs(self, rvecs):
        # Find cell vectors in LAMMPS format
        give_lower(rvecs, self.rvecs)
        self.cell.update_rvecs(self.rvecs)
        # Update the corresponding rotation matrix
        self.update_rot()
        if self.triclinic:
            # Avoid extremely tilted boxes
            xy = self.rvecs[1,0]
            xz = self.rvecs[2,0]
            yz = self.rvecs[2,1]
            x = self.rvecs[0,0]
            y = self.rvecs[1,1]
            z = self.rvecs[2,2]
            #print x,y,z,xy,xz,yz, np.floor( (xy+0.5*x)/x )
            xy -= np.floor( (xy+0.5*x)/x )*x
            xz -= np.floor( (xz+0.5*x)/x )*x
            yz -= np.floor( (yz+0.5*y)/y )*y
            #print np.floor( (xy+0.5*x)/x )
            #print np.floor( (xz+0.5*x)/x )
            #print np.floor( (yz+0.5*y)/y )

            self.lammps.command("change_box all x final %f %30.20f y final %f %30.20f z final %f %30.20f xy final %30.20f xz final %30.20f yz final %30.20f\n" % 
                (0.0,x,0.0,y, 0.0,z,xy,xz,yz))
            #self.lammps.command("change_box all  xy final %30.20f xz final %30.20f yz final %30.20f x final %f %30.20f y final %f %30.20f z final %f %30.20f\n" % 
            #    (xy,xz,yz,0.0,x,0.0,y, 0.0,z))
        else:
            self.lammps.command("change_box all x final %f %30.20f y final %f %30.20f z final %f %30.20f\n" % 
                (0.0,self.rvecs[0,0],0.0,self.rvecs[1,1], 0.0, self.rvecs[2,2]))

    def _internal_compute(self, gpos, vtens):
        with timer.section("LAMMPS overhead"):
            self.update_rvecs(self.system.cell.rvecs)
            self.update_pos(self.system.pos.copy())
        with timer.section("LAMMPS"):
            self.lammps.command("run 0 post no")
        with timer.section("LAMMPS overhead"):
            energy = self.lammps.extract_variable("eng",None,0)
            if gpos is not None:
                f = self.lammps.gather_atoms("f",1,3)
                buffer = np.core.multiarray.int_asbuffer(
                    ctypes.addressof(f), 8*3*self.system.natom)
                gpos[:] = np.frombuffer(buffer, float).reshape((-1,3))
#                for iatom in xrange(self.system.natom):
#                    for j in xrange(3):
#                        gpos[iatom,j] = f[3*iatom+j]
                gpos[:] = -np.einsum('ij,kj', gpos, self.rot.transpose())
            if vtens is not None:
                w = self.lammps.extract_compute("virial",0,1)
                buffer = np.core.multiarray.int_asbuffer(
                    ctypes.addressof(w.contents), 8*6)
                vtens_lammps = np.frombuffer(buffer, float)
#                vtens_lammps = np.zeros(6)
#                for i in xrange(6):
#                    vtens_lammps[i] = w[i]
                # Lammps gives the virial per volume in pascal, so we have to
                # multiply with some prefactors
                vtens_lammps[:] *= -pascal*self.system.cell.volume
                # The [6x1] vector has to be cast to a symmetric [3x3] tensor
                # Lammps orders the values as [xx,yy,zz,xy,xz,yz]
                vtens[np.triu_indices(3)] = vtens_lammps[[0,3,4,1,5,2]]
                vtens[np.tril_indices(3)] = vtens_lammps[[0,3,1,4,5,2]]
                # Finally we have to compute the effect of the rotation on the
                # the virial tensor to get the values in Yaff coordinates
                vtens[:] = np.dot(self.rot.transpose(),np.dot(vtens[:],self.rot))
        return energy

def write_lammps_table(ff, ff2=None, fn='lammps.table', workdir='.', rmin=0.70*angstrom, rmax=12.0*angstrom, nrows=2500, make_plots=False):
    '''
       Write tables containing intermolecular interactions for lammps.
       A table for every pair of ffatypes is  generated.

       **Arguments:**

       ff
            Yaff ForceField instance

       **Optional arguments:**

       fn
            Filename where tables will be stored

       workdir

            Where to save the generated tables

    '''
    numbers = []
    nffatypes = np.amax(ff.system.ffatype_ids) + 1
    # TODO: check that scaling equals 1 for all pairs
    for i in xrange(nffatypes):
        index0 = np.where(ff.system.ffatype_ids==i)[0][0]
        if np.sum(ff.system.ffatype_ids==i)>1:
            index1 = np.where(ff.system.ffatype_ids==i)[0][1]
            numbers.append( [index0, index1] )
        for j in xrange(i+1,nffatypes):
            index1 = np.where(ff.system.ffatype_ids==j)[0][0]
            numbers.append( [index0, index1] )
    # We only consider one neighbor interaction
    ff.compute()
    ff.nlist.nneigh = 1
    # Construct array of evenly spaced values
    distances = np.linspace(rmin, rmax, nrows)
    ftab = open(os.path.join(workdir,fn),'w')
    ftab.write("# LAMMPS tabulated potential generated by Yaff\n")
    ftab.write("# All quantities in atomic units\n")
    ftab.write("# The names of the tables refer to the ffatype_ids that have to be used in the Yaff system\n")
    ftab.write("#%4s %13s %21s %21s\n" % ("i","d","V","F"))
    # Loop over all atom pairs
    for index0, index1 in numbers:
        energies = []
        for d in distances:
            gposnn = np.zeros(ff.system.pos.shape, float)
            ff.nlist.neighs[0] = (index0, index1, d, 0.0, 0.0, d, 0, 0, 0)
            energy = 0
            for part in ff.parts:
                energy += part.compute(gpos=gposnn)
            row = [d, energy, gposnn[index0,2]]
            energies.append( row )
        energies = np.asarray(energies)
        name = '%03d-%03d' % (ff.system.ffatype_ids[index0],
                              ff.system.ffatype_ids[index1])
        ftab.write("%s\nN %d R %f %f\n\n" % (name, nrows, rmin, rmax))
        for irow, row in enumerate(energies):
            ftab.write("%05d %+13.8f %+21.12f %+21.12f\n" % (irow+1, row[0], row[1], row[2]))
        if make_plots:
            if not os.path.isdir(os.path.join(workdir,'lammps_table_plots')): os.mkdir(os.path.join(workdir,'lammps_table_plots'))
            pt.clf()
            #pt.subplot(2,1,1)
            pt.plot(energies[:,0]/angstrom,energies[:,1]/kjmol)
            pt.yscale('symlog',linthreshy=1.0)
            #pt.gca().set_xticks(np.arange(1.5,9.0,1.5))
            #pt.gca().set_yticks(np.arange(-50,75,12.5))
            #pt.xlim([1.5,7.5])
            #pt.ylim([-50.0,50.0])
            pt.xlabel('d [$\AA$]')
            pt.ylabel('E [kJ/mol]')
            #pt.subplot(2,1,2)
            #pt.plot(energies[:,0]/angstrom,energies[:,2])
            pt.savefig(os.path.join(workdir,'lammps_table_plots','%s.png'%name))

def write_lammps_data(system, fn='lammps.data', triclinic=True):
    '''
        Write information about a Yaff system to a LAMMPS data file
            system
                Yaff system

            fn
                Filename to write the LAMMPS data to

            triclinic
                Boolean, specify whether a triclinic cell will be used during
                the simulation. If the cell is orthogonal, set it to False
                as LAMMPS should run slightly faster.
                Default: True
    '''
    if system.cell.nvec != 3:
        raise ValueError('The system must be 3d periodic for Lammps calculations.')
    if system.ffatypes is None:
        raise ValueError('Atom types need to be defined.')
    if system.bonds is None:
        raise ValueError('Bonds need to be defined')
    if system.charges is None:
        charges = np.zeros((system.natom,))
    else:
        charges = system.charges
    fdat = open(fn,'w')
    fdat.write("Generated by Yaff\n\n%20d atoms\n%20d bonds\n%20d angles \n%20d dihedrals\n%20d impropers\n\n" % (system.natom, system.nbond, 0, 0, 0))
    fdat.write("%20d atom types\n%20d bond types\n%20d angle types\n%20d dihedral types\n%20d improper types\n\n" % (np.amax(system.ffatype_ids) + 1, 1,0,0,0) )
    rvecs = np.zeros((3,3))
    give_lower(system.cell.rvecs, rvecs)
    fdat.write("%30.24f %30.24f xlo xhi\n%30.24f %30.24f ylo yhi\n%30.24f %30.24f zlo zhi\n" % (0.0,rvecs[0,0],0.0,rvecs[1,1],0.0,rvecs[2,2]) )
    if triclinic:
        fdat.write("%30.24f %30.24f %30.24f xy xz yz\n" % (rvecs[1,0],rvecs[2,0],rvecs[2,1]) )
    fdat.write("Atoms\n\n")
    for i in xrange(system.natom):
        fdat.write("%5d %3d %3d %30.24f %30.24f %30.24f %30.24f\n" % (i+1,1,system.ffatype_ids[i]+1, charges[i], system.pos[i,0], system.pos[i,1], system.pos[i,2]) )
    fdat.write("\nBonds\n\n")
    for i in xrange(system.nbond):
        fdat.write("%5d %3d %5d %5d\n" % (i+1,1,system.bonds[i,0]+1, system.bonds[i,1]+1))
    fdat.close()

def give_lower(mat0, mat1):
    '''
        Generate a lower triangular matrix mat1 from mat0, where both matrices
        represent cell vectors. See
        http://lammps.sandia.gov/doc/Section_howto.html#howto_12 for more
        information.
    '''
    A = mat0[0,:]
    B = mat0[1,:]
    C = mat0[2,:]
    mat1[:] = 0.0
    mat1[0,0] = np.linalg.norm(A)
    mat1[1,0] = np.dot(B,A)/mat1[0,0]
    if mat1[1,0] > 0.5*mat1[0,0]: mat1[1,0] -= mat1[0,0]
    mat1[1,1] = np.linalg.norm(np.cross(A,B))/mat1[0,0]
    mat1[2,0] = np.dot(C,A)/mat1[0,0]
    if mat1[2,0] > 0.5*mat1[0,0]: mat1[2,0] -= mat1[0,0]
    mat1[2,1] = (np.dot(B,C) - mat1[1,0]*mat1[2,0])/mat1[1,1]
    if mat1[2,1] > 0.5*mat1[1,1]: mat1[2,1] -= mat1[1,1]
    mat1[2,2] = np.sqrt( np.dot(C,C) - mat1[2,0]**2 - mat1[2,1]**2 )
