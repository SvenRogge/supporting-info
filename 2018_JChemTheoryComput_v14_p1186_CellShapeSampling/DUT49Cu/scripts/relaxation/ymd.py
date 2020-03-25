#!/usr/bin/env python

import numpy as np

from molmod.units import angstrom, kcalmol, kjmol
from yaff import ForcePart,StrainCellDOF,CGOptimizer

class PYDLPOLYPart(ForcePart):
    def __init__(self, system, pd):
        ForcePart.__init__(self, 'part_pydlpoly', system)
        self.system = system
        self.pd = pd

    def update_pos(self, pos):
        self.pd.set_xyz(pos/angstrom)

    def update_rvecs(self, rvecs):
        cell = rvecs/angstrom
        self.pd.set_cell(cell, cell_only=True)

    def _internal_compute(self, gpos, vtens):
        # This is necessary to reset the stress tensor
        dlp.nstep = 1
        # Let pydlpoly know that geometry has changed
        self.update_rvecs(self.system.cell.rvecs.copy())
        self.update_pos(self.system.pos.copy())
        # let pydpoly compute energy and forces
        energy, forces = self.pd.calc_energy_force()
        energy *= kcalmol
        forces *= -kcalmol/angstrom
        if gpos is not None:
            gpos[:] = forces
        if vtens is not None:
            stress = self.pd.get_stress().copy()
            stress *= -kcalmol
            vtens[:] = stress
        return energy

if __name__=='__main__':
    # Setup the pydlpoly force field
    import pydlpoly
    from pydlpoly import dlp
    pd = pydlpoly.pydlpoly("dut49_cp")
    pd.setup(bcond=3,keep=True)
    # Construct a Yaff system
    from molmod.periodic import periodic
    from yaff import ForceField, System
    import h5py as h5
    elements = pd.get_elements()
    numbers = np.array([periodic[element.strip().lower().capitalize()].number for element in elements])
    ffatypes = pd.get_atomtypes()
    pos = pd.get_xyz().copy()*angstrom
    rvecs = pd.get_cell().copy()*angstrom
    system = System(numbers, pos, rvecs=rvecs, ffatypes=ffatypes)
    #system = System.from_file("../DUT49_vol45000.chk")
    # Construct a Yaff force field, wrapping around the pydlpoly force field
    part = PYDLPOLYPart(system, pd)
    ff = ForceField(system, [part])

# Check reference trajectory energies
    """with h5.File('../traj.h5','r') as fh5:
        rvecss = fh5['/trajectory/cell'][:]
        poss = fh5['/trajectory/pos'][:]
        epot = fh5['/trajectory/epot'][:]
        vtenss = fh5['trajectory/vtens'][:]
    nstep = rvecss.shape[0]
    for istep in xrange(nstep):
        rvecs = np.array(rvecss[istep], dtype=float)
        pos = np.array(poss[istep], dtype=float)
        ff.update_rvecs(rvecs)
        ff.update_pos(pos)
        gpos = np.zeros(ff.system.pos.shape)
        vtens = np.zeros((3,3))
        ff.compute(gpos,vtens)

        ff.update_rvecs(rvecs)
        ff.update_pos(pos)
        gpos = np.zeros(ff.system.pos.shape)
        vtens = np.zeros((3,3))
        ff.compute(gpos,vtens)

        print istep
        print "-"*100
        print "%20.12f %20.12f" % (epot[istep]/kjmol,ff.energy/kjmol)
        j = 0
        for i in xrange(3):
            print "%20.12f %20.12f" % (vtens[i,j],vtenss[istep][i,j])
    """



    # Do whatever you want with the Yaff force field
    from yaff.sampling import NHCThermostat, MTKBarostat, TBCombination, VerletScreenLog,  XYZWriter,\
        VerletIntegrator, HDF5Writer,RestartWriter
    from molmod.units import kelvin, pascal, femtosecond
    temp = 300*kelvin
    timestep = 0.5*femtosecond
    press=0e6*pascal
    fh5 = h5.File('traj.h5','w')
    hdf5 = HDF5Writer(fh5, step=100)
    restartf=h5.File('restart.h5',mode='w')
# NPT (P=0)
    hooks = [hdf5]
    dof = StrainCellDOF(ff, gpos_rms=1e-6, dpos_rms=1e-4, grvecs_rms=1e-6, drvecs_rms=1e-4, do_frozen=False)
    GradCon=CGOptimizer(dof,hooks=hooks)
    GradCon.run()
    
    ff.system.to_file("opt.chk")


