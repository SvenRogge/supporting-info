import h5py
import numpy as np
from molmod.units import *
from molmod.constants import *
from yaff import *
import matplotlib.pyplot as plt
import math
#plotsettings()
log.set_level(log.silent)

def plotsettings():
    plt.rc(('xtick','ytick','axes'), labelsize=26.0)
    plt.rcParams['lines.linewidth'] = 2.
    plt.rcParams['axes.linewidth'] = 2.
    plt.rcParams['axes.titlesize'] = 30
    plt.rcParams['xtick.major.size'] = 6
    plt.rcParams['xtick.minor.size'] = 2
    plt.rcParams['ytick.major.size'] = 6
    plt.rcParams['ytick.minor.size'] = 2
    plt.rcParams['axes.labelsize'] = 30
    plt.rcParams['figure.subplot.left']=0.
    plt.rcParams['figure.subplot.bottom']=0.
    plt.rcParams['legend.fontsize'] = 12

def plotsettingsax(ax):
    for line in ax.xaxis.get_ticklines() + ax.yaxis.get_ticklines():
        #line.set_color('green')
        line.set_markersize(5)
        line.set_markeredgewidth(1.5)  
             
def plot2dhistogram(timeseries1,timeseries2,savefig=None,showfig=False,xlabel=None,ylabel=None,barlabel=None):
    plotsettings()
    plt.figure(figsize=(15,10))
    bins = 50
    histogram1,bin1 = np.histogram(timeseries1,bins=np.linspace(-1.5,1.5,num=bins))
    histogram2,bin2 = np.histogram(timeseries2,bins=bins)
    binwidth=bin1[1]-bin1[0]
    binheight=bin2[1]-bin2[0]
    histogram2d , xbin2d, ybin2d = np.histogram2d(timeseries1,timeseries2,bins=(bin1,bin2),normed=True)
    x=np.zeros(len(xbin2d)-1)
    for i in xrange(len(x)):
        x[i]=(xbin2d[i]+xbin2d[i+1])/2.
    y=np.zeros(len(ybin2d)-1)
    for i in xrange(len(y)):
        y[i]=(ybin2d[i]+ybin2d[i+1])/2.
    levels=np.arange(0,7,0.5)
    plt.contourf(y,x,histogram2d*10**3,levels)
    #plt.imshow(histogram2d,interpolation='bilinear',origin='lower')
    plt.axis([y.min(),y.max(),x.min(),x.max()])
    #plt.plot([930,930],[-2.,2.],"r--")
    #plt.plot([1080,1080],[-2.,2.],"r--")
    #plt.plot([1260,1260],[-2.,2.],"r--")
    #plt.plot([1450,1450],[-2.,2.],"r--")
    if barlabel!=None:
        plt.colorbar(label=barlabel,ticks=[0,1,2,3,4,5,6])
    else:
        plt.colorbar(ticks=[0,1,2,3,4,5,6])
    plt.title("2D histogram for 121-supercell MIL53(Al) (updated forcefield 10/2016)")
    plt.xticks(np.array([800,1000,1200,1400])*2.)
    if xlabel!=None:
        plt.xlabel(xlabel)
    if ylabel!=None:
        plt.ylabel(ylabel)    
    if savefig!=None:
        plt.savefig(path+savefig,bbox_inches="tight")
    if showfig:
        plt.show()
    else:
        plt.clf()

    
def readh5(filename):
    file = h5py.File(filename,'r')
    celvecs = np.array(file[u'trajectory'][u'cell'])
    positions = np.array(file[u'trajectory'][u'pos'])
    numbers = np.array(file[u'system'][u'numbers'])
    ffatype_ids = np.array(file[u'system'][u'ffatype_ids'])
    ffatypes=np.array(file[u'system'][u'ffatypes'])
    masses= np.array(file[u'system'][u'masses'])
    bonds=np.array(file[u'system'][u'bonds'])
    volumes = np.array(file[u'trajectory'][u'volume'])
    energies= np.array(file[u'trajectory'][u'epot'])
    print "Original h5 file input"
    print "positions"
    print np.shape(positions)
    print "celvecs"
    print np.shape(celvecs)
    print ""
    return np.array(celvecs),np.array(positions),np.array(volumes),np.array(energies)


#####################################################################################
celvecs=[]
positions=[]
volumes=[]
energies=[]

restart=False

if restart:
    for i in xrange(9):
        path="./run2/"+str(i)+"/"
        celvec,position,volume,energy=readh5(path+"mtd-NPs=0T_rand"+str(i)+"_121.h5")
        if i==0:
            celvecs=celvec
            positions=position
            volumes=volume
            energies=energy
        else:
            celvecs=np.concatenate((celvecs,celvec))
            positions=np.concatenate((positions,position))
            volumes=np.concatenate((volumes,volume))
            energies=np.concatenate((energies,energy))
    celvecs=np.array(celvecs)
    positions=np.array(positions)
    volumes=np.array(volumes)
    energies=np.array(energies)
    
        
    volumerange=[750,1550]
    indexes=[]
    for i in xrange(len(volumes)):
        if volumes[i]<volumerange[0]*2.*angstrom**3 or volumes[i]>volumerange[1]*2.*angstrom**3:
            indexes.append(i)
    print "Lenght of out of volume steps "
    print len(indexes)
    print ""
    celvecs=np.delete(celvecs,indexes,axis=0)
    volumes=np.delete(volumes,indexes,axis=0)
    positions=np.delete(positions,indexes,axis=0)         


    print "Shape celvecs " + str(np.shape(celvecs))
    print "Voorbeeld van een celvec voor rotatie "
    print celvecs[1]
            
    for i in xrange(len(celvecs)):
        U,s,V=np.linalg.svd(celvecs[i])
        S=np.diag(s)
        celvecs[i]=np.dot(celvecs[i],np.dot(U,V).T)
        for j in xrange(len(positions[i])):
            positions[i,j]=np.dot(np.dot(U,V),positions[i,j])
    print "Shape celvecs na rotatie " + str(np.shape(celvecs))
    print "Voorbeeld van een celvec na rotatie "
    print celvecs[1]
    print "Shape positions na rotatie " + str(np.shape(positions))
    celvecs=np.reshape(celvecs,(len(celvecs),3*len(celvecs[0])))

    np.save("celvecs",celvecs)
    np.save("positions",positions)
    np.save("volumes",volumes)
    np.save("energies",energies)

if restart==False:
    celvecs=np.load("celvecs.npy")
    positions=np.load("positions.npy")
    volumes=np.load("volumes.npy")
    energies=np.load("energies.npy")

path=""
for i in [1,2,3,5,6,7]:
  
  plot2dhistogram(celvecs[:,i]/angstrom,volumes/angstrom**3,savefig="fig"+str(i)+"_a_y-vol.png",barlabel="$\hat{p}(V,a_y)$ ($1/10^3 \AA ^4$)",xlabel="$V$ ($\AA^3$)",ylabel="$a_y=b_x$ ($\AA$)")
