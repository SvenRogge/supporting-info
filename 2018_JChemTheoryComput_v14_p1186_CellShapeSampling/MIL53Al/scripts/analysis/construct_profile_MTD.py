
import numpy as np
from numpy.polynomial.polynomial import polyval
from molmod.units import *
from yaff import *
import h5py


def getmeanpressvol(file,startvalue=0):
    pressures = file[u'trajectory'][u'press']
    nppressures = np.array(pressures)
    meanpressure = np.mean(nppressures[startvalue:])
    volume = file[u'trajectory'][u'volume'][0]
    if abs(file[u'trajectory'][u'volume'][0] - file[u'trajectory'][u'volume'][-1])>0.01:
        print "ATTENTION: volumes not the same!!"
    return meanpressure,volume
    
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
    
    
"""
V=[]
P=[]
F=[]
startvalue=200
for i in xrange(1440,3080,10):
    print i
    try:
        file = h5py.File( str(i) + '/md-NVs=0T_121_vol'+str(i)+'.h5','r')    
        press,vol = getmeanpressvol(file,startvalue=200)
        P.append(press)
        V.append(vol)
        file.close()
    except:
        print "Volume not yet available"
print V
print P
"""
V=np.load("V.npy")
P=np.load("P.npy")
  
cP,cov = np.polyfit(V,P,11,cov=True)
P_fit=np.polyval(cP,V)
cF=cP/np.array(np.arange(12,0,-1))
cF=np.concatenate((cF,[0]))
F=-np.polyval(cF,V)        


V_np=V[np.argmin(F)]

P_err=poly_error(V,cov)
F_err=int_error(V,V_np,cov)
        
np.save("V.npy",V)
np.save("P.npy",P)
np.save("P_fit.npy",P_fit)
np.save("F.npy",F) 
np.save("P_err.npy",P_err) 
np.save("F_err.npy",F_err) 

