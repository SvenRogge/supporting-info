#! /usr/bin/env python
#
#   Get optimized structure of lp MIL-53
#



import numpy as np
import h5py
from mpi4py import MPI
import sys
import glob #To use *.chk as filename

#Append lammps folder to path
sys.path.append("ff-lammps") 

from molmod.units import *
from yaff import *

from mylammps import *

# Setup MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nproc = comm.Get_size()
# Set random seed, important to get the same velocities for all processes
np.random.seed(5)
# Turn off logging for all processes, it can be turned on for one selected process later on
log.set_level(log.silent)

# Traditional stuff
chk_list=glob.glob('*.chk')
py_list=glob.glob('*.py')
if len(chk_list)==1:
    system = System.from_file(chk_list[0])
else:
    raise Exception("Zero or more than one chk files, aborting (I don't know which one to pick)")
    
ff = ForceField.generate(system, 'ff-lammps/pars.txt', rcut=15*angstrom, alpha_scale=3.2, gcut_scale=1.5, smooth_ei=True)


#Extra checks
'''ff.compute()
for part in ff.parts:
    print part.name, part.energy
print ff.energy'''

# Replace non-covalent part with LAMMPS
for part in ff.parts:
    if part.name=='valence': part_valence = part
    elif part.name=='pair_ei': part_ei = part
# Write LAMMPS system file
write_lammps_data(system,fn='data-lammps/lammps.data')
nlist = BondedNeighborList(system)
pair_pot = PairPotEI(system.charges, 0.0, part_ei.pair_pot.rcut, tr=Switch3(part_ei.pair_pot.get_truncation().width), radii=system.radii)
scalings = Scalings(system, scale1=1.0, scale2=1.0, scale3=0.0)
pair_gauss = ForcePartPair(system,nlist,scalings,pair_pot)
pair_lammps = ForcePartLammps(system, pppm_accuracy=1e-5,fn_system='data-lammps/lammps.data',fn_table='data-lammps/lammps_smoothei2.table', fn_log='none',scalings=np.array([0.0,0.0,1.0,0.0,0.0,1.0]),comm=comm)
ff = ForceField(system,[part_valence,pair_lammps,pair_gauss],nlist)

#Extra checks
'''ff.compute()
for part in ff.parts:
    print part.name, part.energy
print ff.energy'''

#Traditional stuff

temp = 300*kelvin
timestep = 0.5*femtosecond
press = 1e6*pascal

thermo = NHCThermostat(temp, timecon=100.0*femtosecond)
baro = MTKBarostat(ff, temp, press, timecon=1000.0*femtosecond, vol_constraint=True)
TBC = TBCombination(thermo, baro)

hooks=[TBC]
if rank==0:
    f=h5py.File(py_list[0][0:-3]+"_"+chk_list[0][6:-4]+'.h5',mode='w')
    hdf=HDF5Writer(f,step=100)
    hooks.append(hdf)

md = VerletIntegrator(ff, timestep, hooks=hooks)
if len(sys.argv)>1:
    if rank==0:
        print "Running ", sys.argv[1] ," simulation steps."
    md.run(int(sys.argv[1]))
else:
    md.run(800)
