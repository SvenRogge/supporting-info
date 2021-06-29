#! /usr/bin/env python 

from molmod.units import *
import h5py
from yaff import *

fn_ff = 'pars.txt'
fn_chk = 'system.chk'

fn_xyz = 'system_opt.xyz'
fn_h5 = 'system_opt.h5'

fn_opt = 'system_opt.chk'

system = System.from_file(fn_chk)
#system.detect_bonds()
#sys = system.supercell(1,2,1)
ff = ForceField.generate(system, fn_ff, rcut=12*angstrom, smooth_ei=False, gcut_scale=1.5, alpha_scale=3.2)

#dof = CartesianDOF(ff, gpos_rms=1e-8, dpos_rms=1e-6)
dof = StrainCellDOF(ff)

#setup output 
xyz = XYZWriter(fn_xyz)
h5f = h5py.File(fn_h5, 'w')
h5 = HDF5Writer(h5f)
screen = OptScreenLog(step=1) 
opt = CGOptimizer(dof, hooks=[xyz,h5,screen])

#run optimization 
try:
    opt.run()
except KeyboardInterrupt:
    interrupted = True
    log("INTERRUPTED")
    log.hline()

system.to_file(fn_opt)
