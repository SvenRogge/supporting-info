import h5py
import numpy as np
from molmod.units import *
import matplotlib.pyplot as plt
from yaff import *

def readh5(h5filename):
    file = h5py.File(h5filename,'r')
    return file



#Write out init file for system 
log.set_level(log.silent)
def MakeInit(file,systemfilename,initfilename,startvalue=0 ,scalevolume=False,idealvolume=1):
    system = System.from_file(systemfilename)
    print system.cell.rvecs
    print system.pos 
    #Rotate original system unit cell and atomic positions
    U,s,V=np.linalg.svd(system.cell.rvecs)
    S=np.diag(s)
    R=np.dot(U,V)
    #Rotate unit cell vectors
    rvecs=np.dot(system.cell.rvecs,R.T)
    pos=[]
    for j in xrange(len(system.pos)):
        #Rotate nucleic Cartesion positions
        pos.append(np.dot(R,system.pos[j]))
    pos=np.array(pos)
    print rvecs
    print pos
    
    sign_ay=np.sign(rvecs[0,1])
           
    rvecs_h5=file[u'trajectory'][u'cell'][startvalue:,:,:]
    for i,rvec_h5 in enumerate(rvecs_h5):
        #Rotate original system unit cell and atomic positions
        U,s,V=np.linalg.svd(rvec_h5)
        S=np.diag(s)
        R=np.dot(U,V)
        #Rotate unit cell vectors
        rvec_h5=np.dot(rvec_h5,R.T)
        rvecs_h5[i]=rvec_h5
    
    rvecs_h5[:,0,1]=sign_ay*np.abs(rvecs_h5[:,0,1])  
    rvecs_h5[:,1,0]=sign_ay*np.abs(rvecs_h5[:,1,0])     
    
    rvecs_av= np.mean(rvecs_h5,axis=0)
    
    realvolume=np.linalg.det(rvecs_av)
    print "The volume is " +str(realvolume/angstrom**3) + " Angstrom**3"

    if scalevolume==True:
        scaleparameter = (idealvolume/realvolume)**(1./3.)
        rvecs_av = rvecs_av*scaleparameter
    
    print rvecs_av
    system=System(system.numbers,pos,bonds = system.bonds,ffatype_ids=system.ffatype_ids,ffatypes=system.ffatypes,masses=system.masses,rvecs=rvecs_av)
    system.to_file(initfilename)

for volume in xrange(144,308):
    i = volume
    h5filename = str(i)+'0/md-NVs=0T_121_vol'+str(i)+'0.h5'
    file = readh5(h5filename)
    MakeInit(file,str(i)+'0/MIL53_121_vol'+str(volume)+'0.chk','chk_averaged/MIL53_121_vol'+str(i)+'0_averaged_abs.chk',startvalue=100,scalevolume=True,idealvolume=volume*10*angstrom**3)



