''' Ruben Demuynck + Molden Acknowledgement + Yaff Acknowledgement'''
'''This is a non exhaustive list of collective variables, add your own  ... I did not test all cvs in the list'''

import numpy as np
from yaff import *
from molmod.ic import *

class CollectiveVariable(object):
    def __init__(self,name,type,atoms=None,parameter=None):
        self.name=name
        self.type=type
        self.atoms=atoms
        self.parameter=parameter
    def get_value(self,system):
        raise NotImplementedError
    def get_deriv(self,system):
        raise NotImplementedError
    def get_force(self,g,gpos,vtens,system):
        raise NotImplementedError

class Angle(CollectiveVariable):
    def __init__(self,name,atoms):
        if len(atoms)!=3:
            raise ValueError('Distance is defined by two atoms, Angle three atoms and dihedral four atoms')
        CollectiveVariable.__init__(self,name,'angle',atoms=atoms)
    def get_value(self,system):
        return bend_angle(np.array([system.pos[self.atoms[0]],system.pos[self.atoms[1]],system.pos[self.atoms[2]]]))[0]
    def get_deriv(self,system):
        return bend_angle(np.array([system.pos[self.atoms[0]],system.pos[self.atoms[1]],system.pos[self.atoms[2]]]),deriv=1)[1]
    def get_force(self,g,gpos,vtens,system):
        cvderiv=self.get_deriv(system)
        for i,a in enumerate(self.atoms):
                                gpos[a,:]+=g*cvderiv[i,:]
        return gpos,vtens

class Distance(CollectiveVariable):
    def __init__(self,name,atoms):
        if len(atoms)!=2:
            raise ValueError('Distance is defined by two atoms, Angle three atoms and dihedral four atoms')
        CollectiveVariable.__init__(self,name,'distance',atoms=atoms)
    def get_value(self,system):
        return bond_length(np.array([system.pos[self.atoms[0]],system.pos[self.atoms[1]]]))[0]
    def get_deriv(self,system):
        return bond_length(np.array([system.pos[self.atoms[0]],system.pos[self.atoms[1]]]),deriv=1)[1]
    def get_force(self,g,gpos,vtens,system):
        cvderiv=self.get_deriv(system)
        for i,a in enumerate(self.atoms):
                                gpos[a,:]+=g*cvderiv[i,:]
        return gpos,vtens

class Volume(CollectiveVariable):
    def __init__(self):
        CollectiveVariable.__init__(self,'volume','volume')
    def get_value(self,system):
        return system.cell.volume
    def get_force(self,g,gpos,vtens,system):
        vtens+=np.identity(3)*g*self.get_value(system)
        return gpos,vtens

class CellParameter(CollectiveVariable):
    def __init__(self,name,parameter):
        if len(parameter)!=2:
            raise ValueError('Two indexes for the cell parameters required')
        CollectiveVariable.__init__(self,name,'cell',parameter=parameter)
    def get_value(self,system):
        return system.cell.rvecs[self.parameter[0],self.parameter[1]]
    def get_force(self,g,gpos,vtens,system):
        gvec=np.zeros((3,3))
        gvec[self.parameter[0],self.parameter[1]]=g
        vtens+=np.dot(system.cell.rvecs.T,gvec)
        return gpos,vtens
        
        
class CellParameterSymm(CollectiveVariable):
    def __init__(self,name,parameter):
        if len(parameter)!=2:
                        raise ValueError('Two indexes for the cell parameters required')
        CollectiveVariable.__init__(self,name,'cell',parameter=parameter)
    def get_value(self,system):
        rvecs=system.cell.rvecs       
        U,s,V=np.linalg.svd(system.cell.rvecs)
        h_s=np.dot(rvecs,np.dot(U,V).T)
        return h_s[self.parameter[0],self.parameter[1]]
    def get_force(self,g,gpos,vtens,system):
        rvecs=system.cell.rvecs
        U,s,V=np.linalg.svd(rvecs)
        gvec=np.zeros((3,3))
        gvec[self.parameter[0],self.parameter[1]]=g
        gvec=np.dot(gvec,np.dot(U,V).T)
        vtens+=np.dot(system.cell.rvecs.T,gvec)
        return gpos,vtens
