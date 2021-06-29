#!/usr/bin/env python

import numpy as np
import os
from glob import glob

from molmod import UnitCell
from molmod.units import angstrom, deg
from molmod.periodic import periodic
from yaff import System, log
log.set_level(log.silent)

from tools import load_poscar

def load_pdb(filename):
    """Loads a single molecule from a pdb file.

       This function does support only a small fragment from the pdb specification.
       It assumes that there is only one molecular geometry in the pdb file.
    """
    f = file(filename,'r')
    numbers = []
    coordinates = []
    occupancies = []
    betas = []
    data = []
    for line in f:
        if line.startswith('MODEL'):
            data.append([])
        elif line.startswith("ATOM"):
            symbol = line[70:79].strip()
            if symbol=='M':
                number = 99
            else: number = periodic[symbol].number
            numbers.append(number)
            coordinates.append([float(line[30:38])*angstrom, float(line[38:46])*angstrom, float(line[46:54])*angstrom])
            occupancies.append(float(line[54:62]))
        elif line.startswith("CRYST1"):
            lengths = [float(w)*angstrom for w in line.split()[1:4]]
            angles = [float(w)*deg for w in line.split()[4:7]]
            cell = UnitCell.from_parameters3(lengths, angles)
            rvecs = cell.matrix.T.copy()
        elif line.startswith('ENDMDL'):
            data[-1] = [np.array(numbers), np.asarray(coordinates), rvecs]
            numbers = []
            coordinates = []
    f.close()
    if len(data)==0:
        return numbers, np.asarray(coordinates), rvecs
    else:
        return data

def write_filled(fn):
    fw = fn.split(os.sep)[2]
    nmol = int(fn.split(os.sep)[3].split('_')[1])
    irun = int(fn.split(os.sep)[3].split('_')[0].split('run')[1])
    print fn, fw, nmol
    if not os.path.isfile(fn): return
    data = load_pdb(fn)
    numbers, pos, rvecs = data[0]
    ffatypes = ['CCH4','HCH4','HCH4','HCH4','HCH4']*(len(numbers)/5)
    system = System.from_file('../input/%s.chk'%fw)
    natom = system.natom
    newnumbers = np.concatenate((system.numbers, numbers))
    newpos = np.concatenate((system.pos, pos), axis=0)
    newffas = np.concatenate((system.ffatypes[system.ffatype_ids],ffatypes))
    bonds = list(system.bonds)
    for imol in xrange(nmol):
        for i in xrange(4):
            bonds.append([imol*5+i+1+system.natom,imol*5+system.natom])
    bonds = np.asarray(bonds)
    system = System(newnumbers, newpos, rvecs=system.cell.rvecs, ffatypes=newffas, bonds=bonds)
    for iconf in xrange(len(data)):
        if iconf%10!=9: continue
        system.pos[natom:] = data[iconf][1]    
        system.to_file('../filled2/%s_%04d_run%02d_%04d.chk'%(fw,nmol,irun,iconf+1))

def write_snapshot(fw,guest,ff='mm3-mbis',t=0,p=0,suffix='',overwrite=False):
    wdir = os.path.join(fw,ff,guest,'t%d-p%d%s'%(t,p,suffix))
    pattern = os.path.join(wdir,'Movies','System_0','Movie*_component_%s_0.pdb'%guest)
    print pattern
    fns = glob(pattern)
    if len(fns)==0: return
    assert len(fns)==1
    # Check whether data has been added since last time
    if not overwrite:
        snaps = sorted(glob( os.path.join(wdir,'snapshot_*.chk')))
        mtime = max([0] + [os.path.getmtime(snap) for snap in snaps])
        if os.path.getmtime(fns[0])< mtime:
            return 
    print fns[0]
    data = load_pdb(fns[0])
    # Load the framework
    fn_fw = "../ffpars/%s/%s/optcell.chk"%(ff,fw)
    system_fw = System.from_file(fn_fw)
    # Make supercell RASPA used
    supercell = [int(n) for n in fns[0].split('Movie_%s'%fw)[1].split('_')[1].split('.')]
    system_fw = system_fw.supercell(supercell[0],supercell[1],supercell[2])
    # Update positions from RASPA
    fn_fw = os.path.join(wdir,'Movies','System_0','Framework_0_initial.pdb')
    numbers,coordinates,rvecs = load_pdb(fn_fw)
    system_fw.cell.update_rvecs(rvecs)
    #system_fw = System(np.asarray(numbers, dtype=int), coordinates, rvecs=rvecs)
    assert np.all(numbers==system_fw.numbers)
    system_fw.pos[:] = coordinates
    # Load the guest
    fn_guest = "%s.chk"%guest
    system_guest = System.from_file(fn_guest)
    guestnumbers, guestpos, guestrvecs = data[0]
    print guestnumbers, guestpos
    # Loop over the snapshots
    print wdir
    for idata in xrange(len(data)):
        numbers, pos, rvecs = data[idata]
        nmol = numbers.shape[0]/system_guest.natom
        assert nmol*system_guest.natom==numbers.shape[0]
        if nmol==0: continue
        newnumbers = np.concatenate((system_fw.numbers, numbers))
        newpos = np.concatenate((system_fw.pos, pos), axis=0)
        newffas = np.concatenate((system_fw.ffatypes[system_fw.ffatype_ids], list(system_guest.ffatypes[system_guest.ffatype_ids])*nmol))
        bonds = list(system_fw.bonds)
        for imol in xrange(nmol):
            for bond in system_guest.bonds:
                offset = system_fw.natom+imol*system_guest.natom
                bonds.append([offset+bond[0],offset+bond[1]])
        bonds = np.asarray(bonds)
        system = System(newnumbers, newpos, rvecs=system_fw.cell.rvecs, ffatypes=newffas, bonds=bonds)
        snap_fn = os.path.join(wdir,'snapshot_%06d_nmol%04d.chk'%(idata,nmol))
        print snap_fn
        system.to_file(snap_fn)


if __name__=='__main__':
    fws=['ZIF8_opengate','ZIF8_closedgate']
    guests= ['water']
    for fw in fws:
        for guest in guests:
            for p in [26]:
                write_snapshot(fw,guest,p=p,suffix='',ff='qff-dreiding-mbis-tip4p')
