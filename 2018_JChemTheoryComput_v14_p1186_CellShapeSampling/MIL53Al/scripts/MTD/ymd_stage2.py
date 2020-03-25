#! /usr/bin/env python
#
#   Metadynamics
#



import numpy as np
import h5py as h5
from mpi4py import MPI
import sys
import glob #To use *.chk as filename

from molmod.units import *
from yaff import *

#Append lammps folder to path
sys.path.append("ff-lammps") 
from mylammps import *

sys.path.append("metadynamics") 
import colvar_senne as colvar
import DMTD_senne_restart as DMTD

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
press = 0e6*pascal

#Add thermostat and barostat to 'hooks'
thermo = NHCThermostat(temp, timecon=100.0*femtosecond)
baro = MTKBarostat(ff, temp, press, timecon=1000.0*femtosecond, vol_constraint=False)
TBC = TBCombination(thermo, baro)
hooks=[TBC]

#For one process, add a hdf-writer to 'hooks'
if rank==0:
    f=h5.File(py_list[0][0:-3]+"_"+chk_list[0][6:-4]+'.h5',mode='w')
    hdf=HDF5Writer(f,step=100)
    hooks.append(hdf)

#Initiate collective variable for metadynamics and type of hills
cv=colvar.Volume()
width=[25*angstrom**3]
height=0.5*kjmol
vhill=[DMTD.Hills([cv],width=width,height=height)]

#Set simulation length
mtdsteps=2 #Number of hills to be placed
mdsteps=1200 #Number of steps between hill placement
if len(sys.argv)>2:
    if rank==0:
        print "Running ", sys.argv[1] ," metadynamics steps of ", sys.argv[2]," md steps."
    mtdsteps=int(sys.argv[1])   
    mdsteps=int(sys.argv[2])

#Run
replica=DMTD.Metadynamics1D(ff, timestep,mtdsteps,mdsteps,vhill,hooks=hooks)
replica.runMeta()

#Write hill results to file
if rank==0:
    h_pos = ff.part_MTD.B
    h_height = ff.part_MTD.H*height
    h_width = np.tile(width,(mtdsteps,1))
    
    #TODO: Check if shape is (1,mtdsteps)
        
    grid1 = np.array(np.arange(550,1750,1)*2*angstrom**3)
    F=np.zeros(len(grid1))
    for i,pos_i in enumerate(grid1):
        F[i]=-np.sum(h_height*np.exp(-(pos_i-h_pos[:,0])**2/2./h_width[:,0]**2))
    
    #h5-file to write 
    g = h5.File('mtd_results.h5','w')
    #write hill data
    hills = g.create_group('hills')    
    hills.create_dataset('pos', data=h_pos)
    hills.create_dataset('height', data=h_height)
    hills.create_dataset('width', data=h_width)
    #write free energy data    
    fep = g.create_group('fep')    
    fep.create_dataset('volumes', data=grid1)
    fep.create_dataset('F', data=np.array(F))    

