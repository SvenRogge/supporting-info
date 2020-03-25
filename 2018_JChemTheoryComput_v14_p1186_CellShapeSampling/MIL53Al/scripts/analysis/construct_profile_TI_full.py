
import numpy as np
from numpy.polynomial.polynomial import polyval
from molmod.units import *
from yaff import *
import h5py


def getmeanpressvol(file,startvalue=0):
    pressures = file[u'trajectory'][u'press']
    nppressures = np.array(pressures)
    meanpressure = np.mean(nppressures[startvalue:])
    energies = file[u'trajectory'][u'etot']
    npenergies = np.array(energies)
    meanenergy = np.mean(npenergies[startvalue:])
    volume = file[u'trajectory'][u'volume'][0]
    if abs(file[u'trajectory'][u'volume'][0] - file[u'trajectory'][u'volume'][-1])>0.01:
        print "ATTENTION: volumes not the same!!"
    return meanpressure,volume,meanenergy
    
def poly_error(x,cov_matrix):
    pows=len(cov_matrix)
    sigma=[]
    for xi in x:
        # [ x**(pows-1) , x**(pows-2) , ... , x , 1 ]
        xi_pow=np.array([xi**(pows-1-i) for i in xrange(pows)])
        xi_pow=xi_pow/1000.
        sigma_temp=np.dot(xi_pow.T,np.dot(cov,xi_pow))**0.5
        sigma.append(sigma_temp*1000.)
    return sigma
    
def int_error(x,x_ref,cov_matrix):
    pows=len(cov_matrix)
    sigma=[]
    # [ x**(pows) , ... , x ]
    x_ref_pow=np.array([x_ref**(pows-i) for i in xrange(pows)])
    for xi in x:
        # [ x**(pows) , ... , x ]
        xi_pow=np.array([xi**(pows-i) for i in xrange(pows)])
        delta_xi_pow=(x_ref_pow-xi_pow)/np.array(np.arange(12,0,-1))
        sigma.append(np.dot(delta_xi_pow.T,np.dot(cov,delta_xi_pow))**0.5)
    return sigma
    

V=[]
P=[]
F=[]
E=[]
for i in xrange(11520,24321,320):
    print i
    try:
        file = h5py.File('../'+ str(i) + '/md-NVs=0T_242_vol'+str(i)+'.h5','r')    
        press,vol,etot = getmeanpressvol(file,startvalue=2500)
        P.append(press)
        V.append(vol)
        E.append(etot)
        file.close()
    except:
        print "Volume not yet available"

#V=np.load("V.npy")
#P=np.load("P.npy")
  
cP,cov = np.polyfit(V,P,11,cov=True)
V_fit=xrange(11520,24321,16)
P_fit=np.polyval(cP,V_fit)
cF=cP/np.array(np.arange(12,0,-1))
cF=np.concatenate((cF,[0]))
F=-np.polyval(cF,V_fit)        


"""
P_err=poly_error(V,cov)
F_err=int_error(V,V_np,cov)
"""        
np.save("V.npy",V)
print "V_fit"
print np.shape(V_fit)
np.save("V_fit.npy",V_fit)
np.save("P.npy",P)
print "P_fit"
print np.shape(P_fit) 
np.save("P_fit.npy",P_fit)
np.save("F.npy",F) 
print np.shape(E)
np.save("E.npy",E)
"""
np.save("P_err.npy",P_err) 
np.save("F_err.npy",F_err) 
""" 
