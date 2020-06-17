#!/usr/bin/env python

import numpy as np
import os
from glob import glob

from horton import Cell
from horton.units import angstrom, deg
from horton.periodic import periodic

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
            try:
                symbol = line[70:79].strip()
                if symbol == 'M': numbers.append(99)
                else:    numbers.append(periodic[symbol].number)
                coordinates.append([float(line[30:38])*angstrom, float(line[38:46])*angstrom, float(line[46:54])*angstrom])
                occupancies.append(float(line[54:62]))
            except:
                data = data[:-1]
                break
        elif line.startswith("CRYST1"):
            lengths = [float(w)*angstrom for w in line.split()[1:4]]
            angles = [float(w)*deg for w in line.split()[4:7]]
            cell = Cell.from_parameters(lengths, angles)
            rvecs = cell.rvecs.copy()
        elif line.startswith('ENDMDL'):
            data[-1] = [np.array(numbers), np.asarray(coordinates), rvecs]
            numbers = []
            coordinates = []
    f.close()
    if len(data)==0:
        return numbers, np.asarray(coordinates), rvecs
    else:
        return data

def process_loading(fw,nmol,ff='qff-dreiding-mbis-tip4p',guest='water'):
    wdir = os.path.join('..','raspa',fw,ff,guest,'insertion_%d'%(nmol))
    pattern = os.path.join(wdir,'Movies','System_0','Movie*_component_%s_0.pdb'%guest)
    fns = glob(pattern)
    if len(fns)==0: return
    assert len(fns)==1
    print fns[0]
    data = load_pdb(fns[0])
    print len(data)
    restart_fn = os.path.join(os.path.dirname(fns[0]),'Component_%s_0.pdb'%guest)
    print restart_fn
    if os.path.isfile(restart_fn):
        restart_data = load_pdb(restart_fn)
        print len(restart_data)
        data += restart_data

    print len(data)

    # Define adsorption sites
    if fw.endswith('111'):
        adsites_frac = np.array([ [0.0,0.0,0.0], [0.5,0.5,0.5] ])
    elif fw.endswith('gate'):
        tmp = np.array([[0.0,0.0,0.0],[0.25,0.25,0.25]])
        adsites_frac = []
        for a in xrange(2):
            for b in xrange(2):
                for c in xrange(2):
                    for frac in tmp:
                        adsites_frac.append(frac+[0.5*a,0.5*b,0.5*c])
        adsites_frac = np.asarray(adsites_frac)
    else: raise NotImplementedError
    counters = np.zeros((adsites_frac.shape[0],))
    for idata in xrange(len(data)):
        numbers, pos, rvecs = data[idata]
        cell = Cell(rvecs)
        adsites_cart = np.zeros(adsites_frac.shape)
        for isite in xrange(adsites_frac.shape[0]): adsites_cart[isite] = cell.to_cart(adsites_frac[isite])
        for iatom in np.where(numbers==8)[0]:
            deltas = np.array([pos[iatom]-adsite for adsite in adsites_cart])
            for isite in xrange(adsites_cart.shape[0]): cell.mic(deltas[isite])
            distances = np.linalg.norm(deltas, axis=1)
            idx = np.argpartition(distances, 1)
            assert np.amin(distances)<10*angstrom
            counters[np.argmin(distances)] += 1
        #break
    counters /= np.sum(counters)
    print "%20s %5d molecules %8d steps:"%(fw,nmol,len(data)),
    for c in counters: print "%4.1f%%"%(c*100.0),
    print "\n", 

if __name__=='__main__':
    for nmol in [20,40,60,80][:]:
        for suffix in ['opengate','closedgate'][:]:
            process_loading(fw='ZIF8_%s_111'%suffix,nmol=nmol)
            process_loading(fw='ZIF8_%s'%suffix,nmol=8*nmol)
