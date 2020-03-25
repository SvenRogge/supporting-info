import numpy as np
import matplotlib.pyplot as plt
from molmod.units import *
from molmod.constants import *


#Import the plotsettings
from Plotsettings import plotsettings

#For dummylines in legend with marker and line in one
from matplotlib.lines import Line2D
def create_dummy_line(**kwds):
    return Line2D([], [], **kwds)



    
col1='#d73027'
col2='#f46d43'
col3='#fdae61'
col4='#abd9e9'
col5='#74add1'
col6='#4575b4'

    
col1='#b2182b'
col2='#ef8a62'
col3='#fddbc7'
col4='#d1e5f0'
col5='#67a9cf'
col6='#2166ac'

col_TIfull = '#6baed6'
col_MTD = '#045a8d'
col_ipol1 = '#d7301f'
col_ipol2 = '#7f0000'
col_mean = '#f16913'
col_snap = '#feb24c'

def calc_abcalphabetagamma(rvecs):
    cells=np.zeros((len(rvecs),6))
    
    norms=np.sqrt(np.einsum('ikj,ikj->ik',rvecs[:,:,:],rvecs[:,:,:]))
    cells[:,0:3]=norms
    
    product=np.einsum('ij,ij->i',rvecs[:,1,:],rvecs[:,2,:])
    product*=1./(norms[:,1]*norms[:,2])
    alpha=np.arccos(product)/np.pi*180
    cells[:,3]=alpha
    
    product=np.einsum('ij,ij->i',rvecs[:,0,:],rvecs[:,2,:])
    product*=1./(norms[:,0]*norms[:,2])
    beta=np.arccos(product)/np.pi*180
    cells[:,4]=beta
    
    product=np.einsum('ij,ij->i',rvecs[:,1,:],rvecs[:,0,:])
    product*=1./(norms[:,1]*norms[:,0])
    gamma=np.arccos(product)/np.pi*180
    cells[:,5]=gamma
    
    return cells

def plot(fit_vol,fit,col,data_vol=[],data=[],every=1,crit_points=[],crit_points_npy=0,label="",linetype="-",pointtype=".",markersize=4.,dashes=None):
    
    #Don't plot error bars any longer
    if False:#crit_points!=[] :
        if (4,)==np.shape(crit_points):
            plt.errorbar(   crit_points[0]/angstrom**3,
                            crit_points[2]/kjmol,
                            yerr=crit_points[3]*3./kjmol,                    
                            fmt=" ",
                            ecolor=col)
        else:
            plt.errorbar(   crit_points[0,:2]/angstrom**3,
                            crit_points[2,1:]/kjmol,
                            yerr=crit_points[3,1:]*3./kjmol,                    
                            xerr=crit_points[1,:2]*3./angstrom**3,
                            fmt=" ",
                            ecolor=col)
            plt.errorbar(   crit_points[0,2]/angstrom**3,
                            crit_points_npy,                
                            xerr=crit_points[1,2]*3./angstrom**3,
                            fmt=" ",
                            ecolor=col)
    if data_vol!=[] and data!=[]:
        plt.plot(   data_vol/angstrom**3,
                    data/kjmol,
                    color=col,
                    linestyle="",
                    marker=pointtype,
                    markersize=markersize,
                    markevery=every)
    if dashes==None:
        plt.plot(   fit_vol/angstrom**3,
                    fit/kjmol,
                    color=col,
                    linestyle=linetype,
                    label=label)
    else:
        plt.plot(   fit_vol/angstrom**3,
                    fit/kjmol,
                    color=col,
                    linestyle=linetype,
                    label=label,
                    dashes=dashes)



##########################################################################################################

###### METADYNAMICS #####

#121 

F_meta      = np.load("121/metadynamics/F_metadynamics.npy")
F_fit_meta  = np.load("121/metadynamics/F_fit_metadynamics.npy")
V_meta      = np.load("121/metadynamics/V_metadynamics.npy")
V_fit_meta  = np.load("121/metadynamics/V_fit_metadynamics.npy")
crit_meta   = np.real(np.load("121/metadynamics/crit_points.npy"))

#2x4x2

F_meta_242      = np.load("242/metadynamics/F_metadynamics.npy")
F_fit_meta_242  = np.load("242/metadynamics/F_fit_metadynamics.npy")
V_meta_242      = np.load("242/metadynamics/V_metadynamics.npy")
V_fit_meta_242  = np.load("242/metadynamics/V_fit_metadynamics.npy")
crit_meta_242   = np.real(np.load("242/metadynamics/crit_points.npy"))

##### THERMODYNAMIC INTEGRATION #####

# NVsT
F_fit_TI1    = np.load("121/thermodynamic_integration/NVsT/F_fit.npy")
P_TI1        = np.load("121/thermodynamic_integration/NVsT/P.npy")
P_fit_TI1    = np.load("121/thermodynamic_integration/NVsT/P_fit.npy")
PTENS_TI1    = np.mean(np.load("121/thermodynamic_integration/NVsT/PTENS.npy"),axis=1)
V_TI1        = np.load("121/thermodynamic_integration/NVsT/V.npy")/2.
V_fit_TI1    = np.load("121/thermodynamic_integration/NVsT/V_fit.npy")
E_TI1        = np.load("121/thermodynamic_integration/NVsT/E.npy")/2.
crit_TI1     = np.real(np.load("121/thermodynamic_integration/NVsT/crit_points.npy"))

# NVsT 2x4x2
F_fit_TI1_242    = np.load("242/thermodynamic_integration/NVsT/F_fit.npy")
P_TI1_242        = np.load("242/thermodynamic_integration/NVsT/P.npy")
PTENS_TI1_242        = np.mean(np.load("242/thermodynamic_integration/NVsT/PTENS.npy"),axis=1)
P_fit_TI1_242    = np.load("242/thermodynamic_integration/NVsT/P_fit.npy")
V_TI1_242        = np.load("242/thermodynamic_integration/NVsT/V.npy")/2./8.
V_fit_TI1_242    = np.load("242/thermodynamic_integration/NVsT/V_fit.npy")
E_TI1_242        = np.load("242/thermodynamic_integration/NVsT/E.npy")/2./8.
crit_TI1_242     = np.real(np.load("242/thermodynamic_integration/NVsT/crit_points.npy"))

# NVhT snapshot 

F_fit_TI2    = np.load("121/thermodynamic_integration/NhT_snapshot/F_fit.npy")
P_TI2        = np.load("121/thermodynamic_integration/NhT_snapshot/P.npy")
PTENS_TI2    = np.mean(np.load("121/thermodynamic_integration/NhT_snapshot/PTENS.npy"),axis=1)
P_fit_TI2    = np.load("121/thermodynamic_integration/NhT_snapshot/P_fit.npy")
RVECS_TI2    = np.load("121/thermodynamic_integration/NhT_snapshot/RVECS.npy")
CELLS_TI2 = calc_abcalphabetagamma(RVECS_TI2)
V_TI2        = np.load("121/thermodynamic_integration/NhT_snapshot/V.npy")/2.
V_fit_TI2    = np.load("121/thermodynamic_integration/NhT_snapshot/V_fit.npy")
crit_TI2     = np.real(np.load("121/thermodynamic_integration/NhT_snapshot/crit_points.npy"))

# NVhT rescaling

F_fit_TI3    = np.load("121/thermodynamic_integration/NhT_rescaling/P_a/F_fit.npy")
P_TI3        = np.load("121/thermodynamic_integration/NhT_rescaling/P_a/P.npy")
P_fit_TI3    = np.load("121/thermodynamic_integration/NhT_rescaling/P_a/P_fit.npy")
V_TI3        = np.load("121/thermodynamic_integration/NhT_rescaling/P_a/V.npy")/2.
V_fit_TI3    = np.load("121/thermodynamic_integration/NhT_rescaling/P_a/V_fit.npy")
E_TI3        = np.load("121/thermodynamic_integration/NhT_rescaling/P_a/E.npy")/2.
crit_TI3     = np.real(np.load("121/thermodynamic_integration/NhT_rescaling/P_a/crit_points.npy"))

F_fit_TI3_old    = np.load("121/thermodynamic_integration/NhT_rescaling/P/F_fit.npy")
P_TI3_old        = np.load("121/thermodynamic_integration/NhT_rescaling/P/P.npy")
P_fit_TI3_old    = np.load("121/thermodynamic_integration/NhT_rescaling/P/P_fit.npy")
V_TI3_old        = np.load("121/thermodynamic_integration/NhT_rescaling/P/V.npy")/2.
V_fit_TI3_old    = np.load("121/thermodynamic_integration/NhT_rescaling/P/V_fit.npy")
E_TI3_old        = np.load("121/thermodynamic_integration/NhT_rescaling/P/E.npy")/2.
crit_TI3_old     = np.real(np.load("121/thermodynamic_integration/NhT_rescaling/P/crit_points.npy"))

# NVhT best guess

F_fit_TI4    = np.load("121/thermodynamic_integration/NhT_best/F_fit.npy")
P_TI4        = np.load("121/thermodynamic_integration/NhT_best/P.npy")
P_fit_TI4    = np.load("121/thermodynamic_integration/NhT_best/P_fit.npy")
PTENS_TI4    = np.mean(np.load("121/thermodynamic_integration/NhT_best/PTENS.npy"),axis=1)
RVECS_TI4        = np.load("121/thermodynamic_integration/NhT_best/RVECS.npy")
CELLS_TI4 = calc_abcalphabetagamma(RVECS_TI4)
V_TI4        = np.load("121/thermodynamic_integration/NhT_best/V.npy")/2.
V_fit_TI4    = np.load("121/thermodynamic_integration/NhT_best/V_fit.npy")
E_TI4        = np.load("121/thermodynamic_integration/NhT_best/E.npy")/2.
crit_TI4     = np.real(np.load("121/thermodynamic_integration/NhT_best/crit_points.npy"))

# NVhT interpolation 2
"""
F_fit_TI5    = np.load("121/thermodynamic_integration/NhT_ipol2/F_fit.npy")
P_TI5        = np.load("121/thermodynamic_integration/NhT_ipol2/P.npy")
P_fit_TI5    = np.load("121/thermodynamic_integration/NhT_ipol2/P_fit.npy")
V_TI5        = np.load("121/thermodynamic_integration/NhT_ipol2/V.npy")/2.
V_fit_TI5    = np.load("121/thermodynamic_integration/NhT_ipol2/V_fit.npy")
E_TI5    = np.load("121/thermodynamic_integration/NhT_ipol2/E.npy")/2.
crit_TI5     = np.real(np.load("121/thermodynamic_integration/NhT_ipol2/crit_points.npy"))
"""
# NhT ipol_new 1

F_fit_TI6    = np.load("121/thermodynamic_integration/NhT_ipol_1_300K/F_fit.npy")
P_TI6        = np.load("121/thermodynamic_integration/NhT_ipol_1_300K/P.npy")
P_ATI_TI6    = np.load("121/thermodynamic_integration/NhT_ipol_1_300K/P_ATI.npy")
P_fit_TI6    = np.load("121/thermodynamic_integration/NhT_ipol_1_300K/P_fit.npy")
V_TI6        = np.load("121/thermodynamic_integration/NhT_ipol_1_300K/V.npy")/2.
V_fit_TI6    = np.load("121/thermodynamic_integration/NhT_ipol_1_300K/V_fit.npy")
E_TI6        = np.load("121/thermodynamic_integration/NhT_ipol_1_300K/E.npy")/2.
crit_TI6     = np.load("121/thermodynamic_integration/NhT_ipol_1_300K/crit_points.npy")
PTENS_TI6     = np.mean(np.load("121/thermodynamic_integration/NhT_ipol_1_300K/PTENS.npy"),axis=1)
RVECS_TI6     = np.load("121/thermodynamic_integration/NhT_ipol_1_300K/RVECS.npy")
CELLS_TI6 = calc_abcalphabetagamma(RVECS_TI6)

# NhT ipol_new 2

F_fit_TI7    = np.load("121/thermodynamic_integration/NhT_ipol_2_300K/F_fit.npy")
P_TI7        = np.load("121/thermodynamic_integration/NhT_ipol_2_300K/P.npy")
P_ATI_TI7    = np.load("121/thermodynamic_integration/NhT_ipol_2_300K/P_ATI.npy")
P_fit_TI7    = np.load("121/thermodynamic_integration/NhT_ipol_2_300K/P_fit.npy")
V_TI7        = np.load("121/thermodynamic_integration/NhT_ipol_2_300K/V.npy")/2.
V_fit_TI7    = np.load("121/thermodynamic_integration/NhT_ipol_2_300K/V_fit.npy")
E_TI7        = np.load("121/thermodynamic_integration/NhT_ipol_2_300K/E.npy")/2.
crit_TI7     = np.load("121/thermodynamic_integration/NhT_ipol_2_300K/crit_points.npy")
PTENS_TI7     = np.mean(np.load("121/thermodynamic_integration/NhT_ipol_2_300K/PTENS.npy"),axis=1)
RVECS_TI7     = np.load("121/thermodynamic_integration/NhT_ipol_2_300K/RVECS.npy")
CELLS_TI7 = calc_abcalphabetagamma(RVECS_TI7)

# NhT best guess 2x4x2
F_fit_TI4_242    = np.load("242/thermodynamic_integration/NhT/mean/F_fit.npy")/8.
P_TI4_242        = np.load("242/thermodynamic_integration/NhT/mean/P.npy")
PTENS_TI4_242        = np.mean(np.load("242/thermodynamic_integration/NhT/mean/PTENS.npy"),axis=1)
P_fit_TI4_242    = np.load("242/thermodynamic_integration/NhT/mean/P_fit.npy")
RVECS_TI4_242    = np.load("242/thermodynamic_integration/NhT/mean/RVECS.npy")/2.
CELLS_TI4_242 = calc_abcalphabetagamma(RVECS_TI4_242)
V_TI4_242        = np.load("242/thermodynamic_integration/NhT/mean/V.npy")/2./8.
V_fit_TI4_242    = np.load("242/thermodynamic_integration/NhT/mean/V_fit.npy")/8.
E_TI4_242        = np.load("242/thermodynamic_integration/NhT/mean/E.npy")/2./8.
crit_TI4_242     = np.real(np.load("242/thermodynamic_integration/NhT/mean/crit_points.npy"))

# NhT ipol_new 2x4x2
F_fit_TI5_242    = np.load("242/thermodynamic_integration/NhT/ipol/F_fit.npy")
P_TI5_242        = np.load("242/thermodynamic_integration/NhT/ipol/P.npy")
PTENS_TI5_242        = np.load("242/thermodynamic_integration/NhT/ipol/PTENS_mean.npy")
RVECS_TI5_242        = np.load("242/thermodynamic_integration/NhT/ipol/RVECS.npy")/2.
CELLS_TI5_242 = calc_abcalphabetagamma(RVECS_TI5_242)
P_ATI_TI5_242        = np.load("242/thermodynamic_integration/NhT/ipol/P_ATI.npy")
P_fit_TI5_242    = np.load("242/thermodynamic_integration/NhT/ipol/P_fit.npy")
V_TI5_242        = np.load("242/thermodynamic_integration/NhT/ipol/V.npy")/2./8.
V_fit_TI5_242    = np.load("242/thermodynamic_integration/NhT/ipol/V_fit.npy")
E_TI5_242        = np.load("242/thermodynamic_integration/NhT/ipol/E.npy")/2./8.
crit_TI5_242     = np.real(np.load("242/thermodynamic_integration/NhT/ipol/crit_points.npy"))

"""
print np.shape(crit_TI5_242)
assert (4,)==np.shape(crit_TI5_242)
if (4,)==np.shape(crit_TI5_242):
    print "Yes"
"""


##### NORMAL MODE ANALYSIS #####

# NVhT best guess

F_NMA1      = np.load("121/nma/NhT_best/F_nma_av3.npy")
TS_NMA1     = np.load("121/nma/NhT_best/TS_nma_av3.npy")
V_NMA1      = np.load("121/nma/NhT_best/V_nma_av3.npy")
F_fit_NMA1  = np.load("121/nma/NhT_best/F_fit_nma_av3.npy")
TS_fit_NMA1 = np.load("121/nma/NhT_best/TS_fit_nma_av3.npy")
V_fit_NMA1  = np.load("121/nma/NhT_best/V_fit_nma_av3.npy")
crit_NMA1   = np.real(np.load("121/nma/NhT_best/crit_points_nma_av3.npy"))

# 2x4x2 NVhT best guess

F_NMA1_242      = np.load("242/nma/NhT_best/F_nma.npy")/8.
TS_NMA1_242     = np.load("242/nma/NhT_best/TS_nma.npy")/8.
V_NMA1_242      = np.load("242/nma/NhT_best/V_nma.npy")/8.
F_fit_NMA1_242  = np.load("242/nma/NhT_best/F_fit_nma.npy")/8.
TS_fit_NMA1_242 = np.load("242/nma/NhT_best/TS_fit_nma.npy")/8.
V_fit_NMA1_242  = np.load("242/nma/NhT_best/V_fit_nma.npy")/8.
crit_NMA1_242   = np.real(np.load("242/nma/NhT_best/crit_points_nma.npy"))/8.


##### PRINCIPAL COMPONENT ANALYSIS #####

# NVhT best guess

F_PCA1      = np.load("121/pca/NhT_best/F_pca_av3.npy")
TS_PCA1     = np.load("121/pca/NhT_best/TS_pca_av3.npy")
V_PCA1      = np.load("121/pca/NhT_best/V_pca_av3.npy")
F_fit_PCA1  = np.load("121/pca/NhT_best/F_fit_pca_av3.npy")
TS_fit_PCA1 = np.load("121/pca/NhT_best/TS_fit_pca_av3.npy")
V_fit_PCA1  = np.load("121/pca/NhT_best/V_fit_pca_av3.npy")
crit_PCA1   = np.real(np.load("121/pca/NhT_best/crit_points_pca_av3.npy"))

F_PCA1     += crit_PCA1[2,1] #Set np minimum to zero
F_fit_PCA1 += crit_PCA1[2,1] #Set np minimum to zero

# NVs=0T 

F_PCA2      = np.load("121/pca/NVsT/F_pca_NVs=0T.npy")
TS_PCA2     = np.load("121/pca/NVsT/TS_pca_NVs=0T.npy")
V_PCA2      = np.load("121/pca/NVsT/V_pca_NVs=0T.npy")
F_fit_PCA2  = np.load("121/pca/NVsT/F_fit_pca_NVs=0T.npy")
TS_fit_PCA2 = np.load("121/pca/NVsT/TS_fit_pca_NVs=0T.npy")
V_fit_PCA2  = np.load("121/pca/NVsT/V_fit_pca_NVs=0T.npy")
crit_PCA2   = np.real(np.load("121/pca/NVsT/crit_points_pca_NVs=0T.npy"))

F_PCA2     += crit_PCA2[2,1] #Set np minimum to zero
F_fit_PCA2 += crit_PCA2[2,1] #Set np minimum to zero

####################################################################################################################################

#               DUT49

################################################

### OPEN PORE ###

# NVsT
F_fit_TI1_DUT    = np.load("DUT49/open_pore/NVsT/TI/F_fit.npy")
P_TI1_DUT        = np.load("DUT49/open_pore/NVsT/TI/P.npy")
P_fit_TI1_DUT    = np.load("DUT49/open_pore/NVsT/TI/P_fit.npy")
PTENS_TI1_DUT    = np.mean(np.load("DUT49/open_pore/NVsT/TI/PTENS.npy"),axis=1)
V_TI1_DUT        = np.load("DUT49/open_pore/NVsT/TI/V.npy")
V_fit_TI1_DUT    = np.load("DUT49/open_pore/NVsT/TI/V_fit.npy")
E_TI1_DUT        = np.load("DUT49/open_pore/NVsT/TI/E.npy")
crit_TI1_DUT     = np.real(np.load("DUT49/open_pore/NVsT/TI/crit_points.npy"))

# MTD 
F_MTD1_DUT    = np.load("DUT49/open_pore/NVsT/MTD/F_MTD.npy")
V_MTD1_DUT    = np.load("DUT49/open_pore/NVsT/MTD/V_MTD.npy")
F_fit_MTD1_DUT    = np.load("DUT49/open_pore/NVsT/MTD/F_fit_metadynamics.npy")
V_fit_MTD1_DUT    = np.load("DUT49/open_pore/NVsT/MTD/V_fit_metadynamics.npy")
crit_MTD1_DUT    = np.load("DUT49/open_pore/NVsT/MTD/crit_points.npy")

# NVsT
F_fit_TI2_DUT    = np.load("DUT49/open_pore/NVhT/TI/F_fit.npy")
P_TI2_DUT        = np.load("DUT49/open_pore/NVhT/TI/P.npy")
P_fit_TI2_DUT    = np.load("DUT49/open_pore/NVhT/TI/P_fit.npy")
PTENS_TI2_DUT    = np.mean(np.load("DUT49/open_pore/NVhT/TI/PTENS.npy"),axis=1)
V_TI2_DUT        = np.load("DUT49/open_pore/NVhT/TI/V.npy")
V_fit_TI2_DUT    = np.load("DUT49/open_pore/NVhT/TI/V_fit.npy")
E_TI2_DUT        = np.load("DUT49/open_pore/NVhT/TI/E.npy")
crit_TI2_DUT     = np.real(np.load("DUT49/open_pore/NVhT/TI/crit_points.npy"))

# MTD cubic
F_fit_MTD2_DUT    = np.load("DUT49/open_pore/NVhT/MTD/F_MTD_cubic.npy")
V_fit_MTD2_DUT    = np.load("DUT49/open_pore/NVhT/MTD/V_MTD_cubic.npy")
F_fit_MTD2_DUT    = np.load("DUT49/open_pore/NVhT/MTD/F_fit_metadynamics_cubic.npy")
V_fit_MTD2_DUT    = np.load("DUT49/open_pore/NVhT/MTD/V_fit_metadynamics_cubic.npy")
crit_MTD2_DUT    = np.load("DUT49/open_pore/NVhT/MTD/crit_points_cubic.npy")

### CLOSED PORE ###

# NVsT
F_fit_TI3_DUT    = np.load("DUT49/closed_pore/NVsT/TI/F_fit.npy")
P_TI3_DUT        = np.load("DUT49/closed_pore/NVsT/TI/P.npy")
P_fit_TI3_DUT    = np.load("DUT49/closed_pore/NVsT/TI/P_fit.npy")
PTENS_TI3_DUT    = np.mean(np.load("DUT49/closed_pore/NVsT/TI/PTENS.npy"),axis=1)
V_TI3_DUT        = np.load("DUT49/closed_pore/NVsT/TI/V.npy")
V_fit_TI3_DUT    = np.load("DUT49/closed_pore/NVsT/TI/V_fit.npy")
E_TI3_DUT        = np.load("DUT49/closed_pore/NVsT/TI/E.npy")
crit_TI3_DUT     = np.real(np.load("DUT49/closed_pore/NVsT/TI/crit_points.npy"))


# NVhT
F_fit_TI4_DUT    = np.load("DUT49/closed_pore/NVhT/TI/F_fit.npy")
P_TI4_DUT        = np.load("DUT49/closed_pore/NVhT/TI/P.npy")
P_fit_TI4_DUT    = np.load("DUT49/closed_pore/NVhT/TI/P_fit.npy")
PTENS_TI4_DUT    = np.mean(np.load("DUT49/closed_pore/NVhT/TI/PTENS.npy"),axis=1)
V_TI4_DUT        = np.load("DUT49/closed_pore/NVhT/TI/V.npy")
V_fit_TI4_DUT    = np.load("DUT49/closed_pore/NVhT/TI/V_fit.npy")
E_TI4_DUT        = np.load("DUT49/closed_pore/NVhT/TI/E.npy")
crit_TI4_DUT     = np.real(np.load("DUT49/closed_pore/NVhT/TI/crit_points.npy"))
















####################################################################################################################################

plotsettings()       

plot(   V_fit_meta,
        F_fit_meta,
        col1,
        crit_points=crit_meta,
        crit_points_npy=-100,
        linetype="-",
        #data_vol=V_meta,
        #data=F_meta,
        every=5,
        label="MTD")
        
plot(   V_fit_TI1,
        F_fit_TI1,
        col6,
        crit_points=crit_TI1,
        crit_points_npy=-100,
        linetype="-",
        label="TI$\mathregular{^{full}}$")

      
plt.ylabel("$F$ (kJ/mol)")
plt.ylim(0,40)
plt.yticks([0,10,20,30,40]) 
#plt.title("TI$\mathregular{^{full}}$ gives the only precise estimate of the free energy")        
plt.legend( loc="upper center",
            handlelength=3,
            ncol=3,bbox_to_anchor=(0.5, 1.3))
#plt.tight_layout()
plt.savefig("F_reference.png")
plt.savefig("F_reference.pdf")
#plt.show()    
plt.clf() 


######################################     

plotsettings()       

plot(   V_fit_meta,
        F_fit_meta,
        col1,
        crit_points=crit_meta,
        crit_points_npy=-100,
        linetype="-",
        #data_vol=V_meta,
        #data=F_meta,
        every=5,
        label="MTD")
        
plot(   V_fit_TI1,
        F_fit_TI1,
        col6,
        crit_points=crit_TI1,
        crit_points_npy=-100,
        linetype="-",
        label="TI$\mathregular{^{full}}$")
"""plot(   V_fit_TI2,
        F_fit_TI2,
        col2,
        crit_points=crit_TI2,
        crit_points_npy=-100,
        linetype="-",
        #dashes=(2,3,2,3),
        label="TI$\mathregular{^{snap}}$")"""
plot(   V_fit_TI6,
        F_fit_TI6,
        col5,
        linetype="-",
        #crit_points=crit_TI3,
        crit_points_npy=-100,
        #dashes=(10,3,3,3),
        label="TI$\mathregular{^{ipol}}$")
     
plot(   V_fit_TI4,
        F_fit_TI4+crit_TI4[2,1],
        col3,
        crit_points=crit_TI4,
        crit_points_npy=-100,
        linetype="-",
        label="TI$\mathregular{^{mean}}$")        

plot(   V_fit_TI7,
        F_fit_TI7,
        col4,
        linetype="-",
        #dashes=(3,5,3,5),
        crit_points=crit_TI7,
        crit_points_npy=-100,
        label="TI$\mathregular{^{ipol2}}$")
      
plt.ylabel("$F$ (kJ/mol)")
plt.ylim(-10,50)
plt.yticks([0,10,20,30,40]) 
#plt.title("TI$\mathregular{^{full}}$ gives the only precise estimate of the free energy")        
plt.legend( loc="upper center",
            handlelength=3,
            ncol=3,bbox_to_anchor=(0.5, 1.3))
#plt.tight_layout()
plt.savefig("F.png")
plt.savefig("F.pdf")
#plt.show()    
plt.clf() 


######################################
plotsettings()

plot(   V_fit_meta,
        F_fit_meta,
        col1,
        crit_points=crit_meta,
        crit_points_npy=-100,
        linetype="-",
        #data_vol=V_meta,
        #data=F_meta,
        every=5,
        label="$\mathregular{MTD_{1x2x1}}$")
        
plot(   V_fit_TI1,
        F_fit_TI1,
        col6,
        crit_points=crit_TI1,
        crit_points_npy=-100,
        linetype="-",
        label="$\mathregular{TI^{full}_{1x2x1}}$")
        
plot(   V_fit_TI4,
        F_fit_TI4+crit_TI4[2,1],
        col3,
        crit_points=crit_TI4,
        crit_points_npy=-100,
        linetype="-",
        label="$\mathregular{TI^{mean}_{1x2x1}}$")

plot(   V_fit_TI6,
        F_fit_TI6,
        col5,
        crit_points=crit_TI6,
        crit_points_npy=-100,
        linetype="-",
        label="$\mathregular{TI^{ipol}_{1x2x1}}$")
        
plot(   V_fit_meta_242,
        F_fit_meta_242, 
        col1,
        crit_points=crit_meta_242,
        crit_points_npy=-100,
        linetype="--",
        #data_vol=V_meta_242,
        #data=F_meta_242,
        every=5,
        label="$\mathregular{MTD_{2x4x2}}$")
        
plot(   V_fit_TI1_242,
        F_fit_TI1_242, 
        col6,
        crit_points=crit_TI1_242,
        crit_points_npy=-100,
        linetype="--",
        #data_vol=V_meta_242,
        #data=F_meta_242,
        every=5,
        label="$\mathregular{TI^{full}_{2x4x2}}$")
        
plot(   V_fit_TI4_242,
        F_fit_TI4_242, 
        col3,
        crit_points=crit_TI4_242,
        crit_points_npy=-100,
        linetype="--",
        #data_vol=V_meta_242,
        #data=F_meta_242,
        every=5,
        label="$\mathregular{TI^{mean}_{2x4x2}}$")
        
plot(   V_fit_TI5_242,
        F_fit_TI5_242, 
        col5,
        crit_points=crit_TI5_242,
        crit_points_npy=-100,
        linetype="--",
        #data_vol=V_meta_242,
        #data=F_meta_242,
        every=5,
        label="$\mathregular{TI^{ipol}_{2x4x2}}$")
        
"""plot(   V_fit_NMA1_242,
        F_fit_NMA1_242-np.min(F_fit_NMA1_242),
        "r",
        linetype="--",
        data_vol=V_NMA1_242,
        data=F_NMA1_242-np.min(F_fit_NMA1_242),
        crit_points=crit_NMA1_242,
        crit_points_npy=-100,
        label="2x4x2 NMA$\mathregular{^{mean}}$") 
"""        
plt.ylabel("$F$ (kJ/mol)")
plt.ylim(0,40)
plt.yticks([0,10,20,30,40]) 
#plt.title("1x2x1 vs 2x4x2 unit cell")        
plt.legend( loc="upper center",
            handlelength=2,
            ncol=3,
            borderpad=0.2,
            borderaxespad=0.25,
            bbox_to_anchor=(0.5, 1.3))
#plt.tight_layout()
plt.savefig("F_combined.png")
plt.savefig("F_combined.pdf")
#plt.show()    
plt.clf()

####################################################################################

######################################
plotsettings()
        
plot(   V_fit_meta_242,
        F_fit_meta_242, 
        col1,
        crit_points=crit_meta_242,
        crit_points_npy=-100,
        linetype="--",
        #data_vol=V_meta_242,
        #data=F_meta_242,
        every=5,
        label="$\mathregular{MTD_{2x4x2}}$")
        
plot(   V_fit_TI1_242,
        F_fit_TI1_242, 
        col6,
        crit_points=crit_TI1_242,
        crit_points_npy=-100,
        linetype="--",
        #data_vol=V_meta_242,
        #data=F_meta_242,
        every=5,
        label="$\mathregular{TI^{full}_{2x4x2}}$")
        
plot(   V_fit_TI4_242,
        F_fit_TI4_242, 
        col3,
        crit_points=crit_TI4_242,
        crit_points_npy=-100,
        linetype="--",
        #data_vol=V_meta_242,
        #data=F_meta_242,
        every=5,
        label="$\mathregular{TI^{mean}_{2x4x2}}$")
        
plot(   V_fit_TI5_242,
        F_fit_TI5_242, 
        col5,
        crit_points=crit_TI5_242,
        crit_points_npy=-100,
        linetype="--",
        #data_vol=V_meta_242,
        #data=F_meta_242,
        every=5,
        label="$\mathregular{TI^{ipol}_{2x4x2}}$")
        
"""plot(   V_fit_NMA1_242,
        F_fit_NMA1_242-np.min(F_fit_NMA1_242),
        "r",
        linetype="--",
        data_vol=V_NMA1_242,
        data=F_NMA1_242-np.min(F_fit_NMA1_242),
        crit_points=crit_NMA1_242,
        crit_points_npy=-100,
        label="2x4x2 NMA$\mathregular{^{mean}}$") 
"""        
plt.ylabel("$F$ (kJ/mol)")
plt.ylim(0,40)
plt.yticks([0,10,20,30,40]) 
#plt.title("1x2x1 vs 2x4x2 unit cell")        
plt.legend( loc="upper center",
            handlelength=2,
            ncol=2,
            borderpad=0.2,
            borderaxespad=0.25,
            bbox_to_anchor=(0.5, 1.3))
#plt.tight_layout()
plt.savefig("F_242.png")
plt.savefig("F_242.pdf")
#plt.show()    
plt.clf()

####################################################################################

plotsettings()
        
plt.plot(   V_TI1/angstrom**3,
        E_TI1/kjmol,
        color=col6,
        linestyle="",
        marker="o",
        markersize=4.,
        label="TI$\mathregular{^{full}}$")

plt.plot(   V_TI4/angstrom**3,
        E_TI4/kjmol, 
        color=col3,
        linestyle="",
        marker="o",
        markersize=4.,
        label="TI$\mathregular{^{mean}}$")
        
        
plt.plot(   V_TI6/angstrom**3,
        E_TI6/kjmol,
        color=col5,
        linestyle="",
        marker="o",
        markersize=4.,
        #dashes=(10,3,3,3),
        label="TI$\mathregular{^{ipol}}$")
        
plt.plot(   V_TI7/angstrom**3,
        E_TI7/kjmol,
        color=col4,
        linestyle="",
        marker="o",
        markersize=4.,
        #dashes=(3,5,3,5),
        label="TI$\mathregular{^{ipol2}}$")     

        
"""        
plot(   V_NMA1_242,
        (F_NMA1_242+TS_NMA1_242)-np.min(F_NMA1_242+TS_NMA1_242),
        "r",
        linetype="--",
        label="2x4x2 NMA$\mathregular{^{mean}}$")
"""

#Include opt 0K results
#plt.plot([803,1497],[-16521,-16476],"bd")

plt.ylabel("$E$ (kJ/mol)")    
plt.legend( loc="upper center",
            numpoints=1,
            ncol=2,
            bbox_to_anchor=(0.5, 1.325))
#plt.ylim(-2,58)
#plt.yticks([0,10,20,30,40,50])
#plt.tight_layout()
plt.savefig("E.png")
#plt.show()
plt.clf()
####################################################################################

plotsettings()
        
plt.plot(   V_TI1_242/angstrom**3,
        E_TI1_242/kjmol,
        color=col6,
        linestyle="",
        marker="s",
        markersize=4.,
        label="$\mathregular{TI^{full}_{2x4x2}}$")

plt.plot(   V_TI4_242/angstrom**3,
        E_TI4_242/kjmol, 
        color=col3,
        linestyle="",
        marker="s",
        markersize=4.,
        label="$\mathregular{TI^{mean}_{2x4x2}}$")
        
        
plt.plot(   V_TI5_242/angstrom**3,
        E_TI5_242/kjmol,
        color=col5,
        linestyle="",
        marker="s",
        markersize=4.,
        #dashes=(10,3,3,3),
        label="$\mathregular{TI^{ipol}_{2x4x2}}$")
         


plt.ylabel("$E$ (kJ/mol)")    
plt.legend( loc="upper center",
            numpoints=1,
            ncol=2,
            bbox_to_anchor=(0.5, 1.35))
#plt.ylim(-2,58)
#plt.yticks([0,10,20,30,40,50])
#plt.tight_layout()
plt.savefig("E_242.png")
plt.clf()
####################################################################################

plotsettings()
"""
plt.plot(   V_TI1/angstrom**3,
            ((E_TI1-F_fit_TI1)-np.array(E_TI1-F_fit_TI1)[20])/kjmol,
            color=col6,
            linestyle="",
            marker="o",
            label="TI$\mathregular{^{full}}$",
            markersize=4.)  
"""       
indexes = np.searchsorted(V_fit_TI1_242,V_TI1_242)
TS_TI1_242 = (E_TI1_242[:-1]-F_fit_TI1_242[indexes[:-1]])-np.array(E_TI1_242[:-1]-F_fit_TI1_242[indexes[:-1]])[5]        
plt.plot(   V_TI1_242[:-1]/angstrom**3,
            (TS_TI1_242)/kjmol,
            color=col6,
            linestyle="",
            marker="s",
            label="$\mathregular{TI^{full}_{2x4x2}}$",
            markersize=4.)   
            
p_temp = np.polyfit(V_fit_TI5_242,F_fit_TI5_242,7)   
F_temp = np.polyval(p_temp,V_TI5_242)        
TS_TI5_242 = (E_TI5_242-F_temp)-np.min(E_TI5_242-F_temp)
plt.plot(   V_TI5_242/angstrom**3,
            (TS_TI5_242)/kjmol,
            color=col5,
            linestyle="",
            marker="s",
            label="$\mathregular{TI^{ipol}_{2x4x2}}$",
            markersize=4)
            
p_temp = np.polyfit(V_fit_TI4_242,F_fit_TI4_242,7)   
F_temp = np.polyval(p_temp,V_TI4_242)        
TS_TI4_242 = (E_TI4_242-F_temp)-E_TI4_242[5]-F_temp[5]
plt.plot(   V_TI4_242/angstrom**3,
            (TS_TI4_242)/kjmol,
            color=col3,
            linestyle="",
            marker="s",
            label="$\mathregular{TI^{mean}_{2x4x2}}$",
            markersize=4)
         

plt.ylabel("$T \Delta S$ (kJ/mol)")    
plt.legend( loc="upper center",
            numpoints=1,
            ncol=2,
            bbox_to_anchor=(0.5, 1.3))
#plt.ylim(-12,42)
plt.yticks([-10,0,10,20,30,40,50,60])
#plt.tight_layout()
plt.savefig("S_242.png")
plt.clf()

        
####################################################################################
'''
plotsettings()
"""
plot(   V_fit_NMA1,
        F_fit_NMA1,
        "r",
        crit_points=crit_NMA1,
        crit_points_npy=-100,
        dashes=(10,3,3,3),
        label="NMA$\mathregular{^{mean}}$")
"""        
plot(   V_fit_TI1,
        F_fit_TI1,
        "b",
        crit_points=crit_TI1,
        crit_points_npy=-100,
        linetype="-",
        label="TI$\mathregular{^{full}}$")

plot(   V_fit_meta,
        F_fit_meta,
        "g",
        crit_points=crit_meta,
        crit_points_npy=-100,
        linetype="-",
        #data_vol=V_meta,
        #data=F_meta,
        every=5,
        label="MTD")


plt.ylabel("$F$ (kJ/mol)")
plt.ylim(-5,35)
plt.yticks([0,10,20,30]) 
#plt.title("Reference free energy profiles")        
plt.legend( loc="upper center",
            handlelength=3,
            ncol=2)
plt.savefig("F_all.png")
#plt.show()    
plt.clf()

'''

########################################################################"

def find_zeros(V,P):
    vol_switch=[]
    for i in xrange(len(P)-1):
        if P[i]*P[i+1]<0:
            vol_switch.append(V[i+1])
    if len(vol_switch)!=0:
        return np.array(vol_switch)/angstrom**3
    else:
        return np.array([0,2000])*angstrom*3
            

def plot_fill_P(ax,V,P,label="",lt="-",col="",alpha=1):
    Vs = find_zeros(V,P)
    plot(V,P*kjmol/pascal/10**6,col,linetype=lt,label=label)
    index=np.logical_and(P>0,V/angstrom**3>np.min(Vs))
    ax.fill_between(V/angstrom**3,P/pascal/10**6,where=index,alpha=alpha,color=col,lw=0)  
    index=np.logical_and(P<0,V/angstrom**3<np.max(Vs))
    ax.fill_between(V/angstrom**3,P/pascal/10**6,where=index,alpha=alpha,color=col,lw=0)  



f,ax=plt.subplots()

plotsettings()
plot_fill_P(ax,V_fit_TI1,P_fit_TI1,label="TI$\mathregular{^{full}}$",lt="-",col=col6,alpha=0.5)

plot_fill_P(ax,V_fit_TI4,P_fit_TI4,label="TI$\mathregular{^{mean}}$",lt="-",col=col3,alpha=0.5)


plt.plot(V_TI1/angstrom**3,P_TI1/pascal/10**6,"o",linestyle="",markerfacecolor=col6,markeredgecolor="k",markersize=4.,alpha=0.5)

plt.plot(V_TI4/angstrom**3,P_TI4/pascal/10**6,"o",linestyle="",markerfacecolor=col3,markeredgecolor="k",markersize=4.,alpha=0.5)


#plot_fill_P(ax,V_fit_TI1_242,P_fit_TI1_242,label="$\mathregular{TI^{full}_{2x4x2}}$",lt="--",col="b",alpha=0.25)

#plot_fill_P(ax,V_fit_TI4_242,P_fit_TI4_242,label="$\mathregular{TI^{mean}_{2x4x2}}$",lt="--",col="m",alpha=0.4)

plt.ylabel("$P$ (MPa)")
plt.ylim(-0.8e3,0.8e3)
plt.yticks([-0.5e3,0,0.5e3]) 
#plt.title("1x2x1 vs 2x4x2 unit cell")
plt.xlabel("$V$ ($\AA ^3$)")
plt.xticks([750,1000,1250,1500])
plt.xlim(700,1550)
plt.legend(loc="upper center",ncol=2,numpoints=1,bbox_to_anchor=(0.5, 1.3))
#plt.tight_layout()
plt.savefig("P_2.png")
plt.savefig("P_2.pdf")
#plt.show()
plt.clf()

########################################

plotsettings()

f,ax=plt.subplots()

plotsettings()


plot_fill_P(ax,V_fit_TI1_242,P_fit_TI1_242,label="$\mathregular{TI^{full}_{2x4x2}}$",lt="--",col=col6,alpha=0.5)

plot_fill_P(ax,V_fit_TI4_242,P_fit_TI4_242,label="$\mathregular{TI^{mean}_{2x4x2}}$",lt="--",col=col3,alpha=0.5)

plt.plot(V_TI1_242/angstrom**3,P_TI1_242/pascal/10**6,"s",linestyle="",markerfacecolor=col6,markeredgecolor="k",markersize=4.,alpha=1)

plt.plot(V_TI4_242/angstrom**3,P_TI4_242/pascal/10**6,"s",linestyle="",markerfacecolor=col3,markeredgecolor="k",markersize=4.,alpha=1)


plt.ylabel("$P$ (MPa)")
plt.ylim(-0.8e3,0.8e3)
plt.yticks([-0.5e3,0,0.5e3]) 
#plt.title("1x2x1 vs 2x4x2 unit cell")
plt.xlabel("$V$ ($\AA ^3$)")
plt.xticks([750,1000,1250,1500])
plt.xlim(700,1550)
plt.legend(loc="upper center",ncol=2,numpoints=1,bbox_to_anchor=(0.5, 1.3))
#plt.tight_layout()
plt.savefig("P_3.png")
plt.savefig("P_3.pdf")
plt.savefig("P_3.svg")
#plt.show()
plt.clf()

########################################
plotsettings()

f,ax=plt.subplots()




p_temp = np.polyfit(V_TI6,P_TI6,7)
plot_fill_P(ax,V_TI6,np.polyval(p_temp,V_TI6),lt="-",col=col5,alpha=0.)
plot_fill_P(ax,V_fit_TI6,P_fit_TI6,lt="-",col=col5,alpha=1.)


plt.plot(V_TI6/angstrom**3,P_TI6/pascal/10**6,color=col5,marker="o",linestyle="",markersize=4.,label="$P$")
plt.plot(V_TI6/angstrom**3,P_ATI_TI6/pascal/10**6,color=col5,marker="o",markerfacecolor="None",alpha=1.,linestyle="",markersize=4.,label="$P_a$")
plt.plot(V_TI6/angstrom**3,P_ATI_TI6/pascal/10**6,color=col5,marker="o",markerfacecolor="white",alpha=0.5,linestyle="",markersize=4.)

plt.title("TI$\mathregular{^{ipol}}$")
plt.ylabel("$P$ (MPa)")
plt.ylim(-0.8e3,0.8e3)
plt.yticks([-0.5e3,0,0.5e3]) 
#plt.title("1x2x1 vs 2x4x2 unit cell")
plt.xlabel("$V$ ($\AA ^3$)")
plt.xticks([750,1000,1250,1500])
plt.xlim(700,1550)
plt.legend(loc="upper center",ncol=2,numpoints=1,bbox_to_anchor=(0.5, 1.3))
#plt.tight_layout()
plt.savefig("P_4.png")
plt.savefig("P_4.pdf")
#plt.show()
plt.clf()

########################################
plotsettings()

f,ax=plt.subplots()

p_temp = np.polyfit(V_TI7,P_TI7,7)
plot_fill_P(ax,V_fit_TI7,np.polyval(p_temp,V_fit_TI7),lt="-",col=col4,alpha=0.)
plot_fill_P(ax,V_fit_TI7,P_fit_TI7,lt="-",col=col4,alpha=1)

plt.plot(V_TI7/angstrom**3,P_TI7/pascal/10**6,color=col4,marker="o",linestyle="",markersize=4.,alpha=1,label="$P$")
plt.plot(V_TI7/angstrom**3,P_ATI_TI7/pascal/10**6,color=col4,marker="o",markerfacecolor="None",alpha=1.,linestyle="",markersize=4.,label="$P_a$")
plt.plot(V_TI7/angstrom**3,P_ATI_TI7/pascal/10**6,color=col4,marker="o",markerfacecolor="white",alpha=0.5,linestyle="",markersize=4.)

plt.title("TI$\mathregular{^{ipol2}}$")
plt.ylabel("$P$ (MPa)")
plt.ylim(-0.8e3,0.8e3)
plt.yticks([-0.5e3,0,0.5e3]) 
#plt.title("1x2x1 vs 2x4x2 unit cell")
plt.xlabel("$V$ ($\AA ^3$)")
plt.xticks([750,1000,1250,1500])
plt.xlim(700,1550)
plt.legend(loc="upper center",ncol=2,numpoints=1,bbox_to_anchor=(0.5, 1.3))
#plt.tight_layout()
plt.savefig("P_5.png")
plt.savefig("P_5.pdf")
#plt.show()
plt.clf()

########################################
plotsettings()

f,ax=plt.subplots()

p_temp = np.polyfit(V_TI5_242,P_ATI_TI5_242,7)
plot_fill_P(ax,V_fit_TI5_242,np.polyval(p_temp,V_fit_TI5_242),lt="-",col=col5,alpha=0.)
plot_fill_P(ax,V_fit_TI5_242,np.polyval(p_temp,V_fit_TI5_242),lt="-",col=col5,alpha=1)

plt.plot(V_TI5_242/angstrom**3,P_TI5_242/pascal/10**6,color=col5,marker="s",linestyle="",markersize=4.,alpha=1,label="$P$")
plt.plot(V_TI5_242/angstrom**3,P_ATI_TI5_242/pascal/10**6,color=col5,marker="s",markerfacecolor="None",alpha=1.,linestyle="",markersize=4.,label="$P_a$")
plt.plot(V_TI5_242/angstrom**3,P_ATI_TI5_242/pascal/10**6,color=col4,marker="s",markerfacecolor="white",alpha=0.5,linestyle="",markersize=4.)

plt.title("$\mathregular{TI^{ipol}_{2x4x2}}$")
plt.ylabel("$P$ (MPa)")
plt.ylim(-0.8e3,0.8e3)
plt.yticks([-0.5e3,0,0.5e3]) 
#plt.title("1x2x1 vs 2x4x2 unit cell")
plt.xlabel("$V$ ($\AA ^3$)")
plt.xticks([750,1000,1250,1500])
plt.xlim(700,1550)
plt.legend(loc="upper center",ncol=2,numpoints=1,bbox_to_anchor=(0.5, 1.15))
#plt.tight_layout()
plt.savefig("P_6.png")
plt.savefig("P_6.pdf")
#plt.show()
plt.clf()

####################################################################################

f,ax=plt.subplots()

plotsettings()
plot_fill_P(ax,V_fit_TI1,P_fit_TI1,label="TI$\mathregular{^{full}}$",lt="-",col=col6,alpha=0.5)

plot_fill_P(ax,V_fit_TI2,P_fit_TI2,label="TI$\mathregular{^{snap}}$",lt="-",col=col2,alpha=0.5)

plt.plot(V_TI1/angstrom**3,P_TI1/pascal/10**6,"o",linestyle="",markerfacecolor=col6,markeredgecolor="k",markersize=4.,alpha=0.5)

plt.plot(V_TI2/angstrom**3,P_TI2/pascal/10**6,"o",linestyle="",markerfacecolor=col2,markeredgecolor="k",markersize=4.,alpha=0.5)

#plot_fill_P(ax,V_fit_TI1_242,P_fit_TI1_242,label="$\mathregular{TI^{full}_{2x4x2}}$",lt="--",col="b",alpha=0.25)

#plot_fill_P(ax,V_fit_TI4_242,P_fit_TI4_242,label="$\mathregular{TI^{mean}_{2x4x2}}$",lt="--",col="m",alpha=0.4)

plt.ylabel("$P$ (MPa)")
plt.ylim(-0.8e3,0.8e3)
plt.yticks([-0.5e3,0,0.5e3]) 
#plt.title("1x2x1 vs 2x4x2 unit cell")
plt.xlabel("$V$ ($\AA ^3$)")
plt.xticks([750,1000,1250,1500])
plt.xlim(700,1550)
plt.legend(loc="upper center",ncol=2,numpoints=1,bbox_to_anchor=(0.5, 1.3))
#plt.tight_layout()
plt.savefig("P_7.png")
plt.savefig("P_7.pdf")
#plt.show()
plt.clf()
        
####################################################################################

plotsettings()

linestyles=["-","-.","--"]
markers=["+","x","1"]
for i in xrange(3):
    for j in xrange(3):
        if i == j: 
            plt.plot(V_TI1/angstrom**3,PTENS_TI1[:,i,j]/pascal/10**6,color=col6,lw=1.5,linestyle=linestyles[i],label=str(i)+str(j))
        else:
            plt.plot(V_TI1/angstrom**3,PTENS_TI1[:,i,j]/pascal/10**6,color="k",lw=1.5,linestyle=linestyles[i+j-1],alpha=1,label=str(i)+str(j))


plt.title("TI$\mathregular{^{full}}$")
plt.ylabel("$\sigma$ (MPa)")
plt.ylim(-1.5e3,1.5e3)
plt.yticks([-1e3,-0.5e3,0,0.5e3,1e3]) 
#plt.title("1x2x1 vs 2x4x2 unit cell")
plt.xlabel("$V$ ($\AA ^3$)")
plt.xticks([750,1000,1250,1500])
plt.xlim(700,1550)
plt.legend(loc="upper center",ncol=3,numpoints=1,bbox_to_anchor=(0.5, 1.15))
#plt.tight_layout()
plt.savefig("stress/full.png")
plt.savefig("stress/full.pdf")
#plt.show()
plt.clf()

####################################################################################

plotsettings()

linestyles=["-","-.","--"]
markers=["+","x","1"]
for i in xrange(3):
    for j in xrange(3):
        if i == j: 
            plt.plot(V_TI4/angstrom**3,PTENS_TI4[:,i,j]/pascal/10**6,color=col3,lw=1.5,linestyle=linestyles[i],label=str(i)+str(j))
        else:
            plt.plot(V_TI4/angstrom**3,PTENS_TI4[:,i,j]/pascal/10**6,color="k",lw=1.5,linestyle=linestyles[i+j-1],alpha=1,label=str(i)+str(j))


plt.title("TI$\mathregular{^{mean}}$ 300K")
plt.ylabel("$\sigma$ (MPa)")
plt.ylim(-1.5e3,1.5e3)
plt.yticks([-1e3,-0.5e3,0,0.5e3,1e3]) 
#plt.title("1x2x1 vs 2x4x2 unit cell")
plt.xlabel("$V$ ($\AA ^3$)")
plt.xticks([750,1000,1250,1500])
plt.xlim(700,1550)
plt.legend(loc="upper center",ncol=3,numpoints=1,bbox_to_anchor=(0.5, 1.15))
#plt.tight_layout()
plt.savefig("stress/mean.png")
plt.savefig("stress/mean.pdf")
#plt.show()
plt.clf()

####################################################################################

plotsettings()

linestyles=["-","-.","--"]
markers=["+","x","1"]
for i in xrange(3):
    for j in xrange(3):
        if i == j: 
            plt.plot(V_TI2/angstrom**3,PTENS_TI2[:,i,j]/pascal/10**6,color=col2,lw=1.5,linestyle=linestyles[i],label=str(i)+str(j))
        else:
            plt.plot(V_TI2/angstrom**3,PTENS_TI2[:,i,j]/pascal/10**6,color="k",lw=1.5,linestyle=linestyles[i+j-1],alpha=1,label=str(i)+str(j))


plt.title("TI$\mathregular{^{snap}}$ 300K")
plt.ylabel("$\sigma$ (MPa)")
plt.ylim(-1.5e3,1.5e3)
plt.yticks([-1e3,-0.5e3,0,0.5e3,1e3]) 
#plt.title("1x2x1 vs 2x4x2 unit cell")
plt.xlabel("$V$ ($\AA ^3$)")
plt.xticks([750,1000,1250,1500])
plt.xlim(700,1550)
plt.legend(loc="upper center",ncol=3,numpoints=1,bbox_to_anchor=(0.5, 1.15))
#plt.tight_layout()
plt.savefig("stress/snap.png")
plt.savefig("stress/snap.pdf")
#plt.show()
plt.clf()
        
####################################################################################

plotsettings()

linestyles=["-","-.","--"]
markers=["+","x","1"]
for i in xrange(3):
    for j in xrange(3):
        if i == j: 
            plt.plot(V_TI6/angstrom**3,PTENS_TI6[:,i,j]/pascal/10**6,color=col5,lw=1.5,linestyle=linestyles[i],label=str(i)+str(j))
        else:
            plt.plot(V_TI6/angstrom**3,PTENS_TI6[:,i,j]/pascal/10**6,color="k",lw=1.5,linestyle=linestyles[i+j-1],alpha=1,label=str(i)+str(j))

plt.title("TI$\mathregular{^{ipol}}$ 300K")
plt.ylabel("$\sigma$ (MPa)")
plt.ylim(-1.5e3,1.5e3)
plt.yticks([-1e3,-0.5e3,0,0.5e3,1e3]) 
#plt.title("1x2x1 vs 2x4x2 unit cell")
plt.xlabel("$V$ ($\AA ^3$)")
plt.xticks([750,1000,1250,1500])
plt.xlim(700,1550)
plt.legend(loc="upper center",ncol=3,numpoints=1,bbox_to_anchor=(0.5, 1.15))
#plt.tight_layout()
plt.savefig("stress/ipol.png")
plt.savefig("stress/ipol.pdf")
#plt.show()
plt.clf()
        
####################################################################################

        

plotsettings()

linestyles=["-","-.","--"]
markers=["+","x","1"]
for i in xrange(3):
    for j in xrange(3):
        if i == j: 
            plt.plot(V_TI7/angstrom**3,PTENS_TI7[:,i,j]/pascal/10**6,color=col4,lw=1.5,linestyle=linestyles[i],label=str(i)+str(j))
        else:
            plt.plot(V_TI7/angstrom**3,PTENS_TI7[:,i,j]/pascal/10**6,color="k",lw=1.5,linestyle=linestyles[i+j-1],alpha=1,label=str(i)+str(j))


plt.title("TI$\mathregular{^{ipol2}}$ 300K")
plt.ylabel("$\sigma$ (MPa)")
plt.ylim(-1.5e3,1.5e3)
plt.yticks([-1e3,-0.5e3,0,0.5e3,1e3]) 
#plt.title("1x2x1 vs 2x4x2 unit cell")
plt.xlabel("$V$ ($\AA ^3$)")
plt.xticks([750,1000,1250,1500])
plt.xlim(700,1550)
plt.legend(loc="upper center",ncol=3,numpoints=1,bbox_to_anchor=(0.5, 1.15))
#plt.tight_layout()
plt.savefig("stress/ipol2.png")
plt.savefig("stress/ipol2.pdf")
#plt.show()
plt.clf()

####################################################################################

        

plotsettings()

linestyles=["-","-.","--"]
markers=["+","x","1"]
for i in xrange(3):
    for j in xrange(3):
        if i == j: 
            plt.plot(V_TI1_242/angstrom**3,PTENS_TI1_242[:,i,j]/pascal/10**6,color=col6,lw=1.5,linestyle=linestyles[i],label=str(i)+str(j))
        else:
            plt.plot(V_TI1_242/angstrom**3,PTENS_TI1_242[:,i,j]/pascal/10**6,color="k",lw=1.5,linestyle=linestyles[i+j-1],alpha=1,label=str(i)+str(j))


plt.title("$\mathregular{TI^{full}_{2x4x2}}$")
plt.ylabel("$\sigma$ (MPa)")
plt.ylim(-1.5e3,1.5e3)
plt.yticks([-1e3,-0.5e3,0,0.5e3,1e3]) 
#plt.title("1x2x1 vs 2x4x2 unit cell")
plt.xlabel("$V$ ($\AA ^3$)")
plt.xticks([750,1000,1250,1500])
plt.xlim(700,1550)
plt.legend(loc="upper center",ncol=3,numpoints=1,bbox_to_anchor=(0.5, 1.15))
#plt.tight_layout()
plt.savefig("stress/full_242.png")
plt.savefig("stress/full_242.pdf")
#plt.show()
plt.clf()

####################################################################################

plotsettings()

linestyles=["-","-.","--"]
markers=["+","x","1"]
for i in xrange(3):
    for j in xrange(3):
        if i == j: 
            plt.plot(V_TI4_242/angstrom**3,PTENS_TI4_242[:,i,j]/pascal/10**6,color=col3,lw=1.5,linestyle=linestyles[i],label=str(i)+str(j))
        else:
            plt.plot(V_TI4_242/angstrom**3,PTENS_TI4_242[:,i,j]/pascal/10**6,color="k",lw=1.5,linestyle=linestyles[i+j-1],alpha=1,label=str(i)+str(j))


plt.title("$\mathregular{TI^{mean}_{2x4x2}}$")
plt.ylabel("$\sigma$ (MPa)")
plt.ylim(-1.5e3,1.5e3)
plt.yticks([-1e3,-0.5e3,0,0.5e3,1e3]) 
#plt.title("1x2x1 vs 2x4x2 unit cell")
plt.xlabel("$V$ ($\AA ^3$)")
plt.xticks([750,1000,1250,1500])
plt.xlim(700,1550)
plt.legend(loc="upper center",ncol=3,numpoints=1,bbox_to_anchor=(0.5, 1.15))
#plt.tight_layout()
plt.savefig("stress/mean_242.png")
plt.savefig("stress/mean_242.pdf")
#plt.show()
plt.clf()

####################################################################################

plotsettings()

linestyles=["-","-.","--"]
markers=["+","x","1"]
for i in xrange(3):
    for j in xrange(3):
        if i == j: 
            plt.plot(V_TI5_242/angstrom**3,PTENS_TI5_242[:,i,j]/pascal/10**6,color=col5,lw=1.5,linestyle=linestyles[i],label=str(i)+str(j))
        else:
            plt.plot(V_TI5_242/angstrom**3,PTENS_TI5_242[:,i,j]/pascal/10**6,color="k",lw=1.5,linestyle=linestyles[i+j-1],alpha=1,label=str(i)+str(j))


plt.title("$\mathregular{TI^{ipol}_{2x4x2}}$")
plt.ylabel("$\sigma$ (MPa)")
plt.ylim(-1.5e3,1.5e3)
plt.yticks([-1e3,-0.5e3,0,0.5e3,1e3]) 
#plt.title("1x2x1 vs 2x4x2 unit cell")
plt.xlabel("$V$ ($\AA ^3$)")
plt.xticks([750,1000,1250,1500])
plt.xlim(700,1550)
plt.legend(loc="upper center",ncol=3,numpoints=1,bbox_to_anchor=(0.5, 1.15))
#plt.tight_layout()
plt.savefig("stress/ipol_242.png")
plt.savefig("stress/ipol_242.pdf")
#plt.show()
plt.clf()

####################################################################################

plotsettings()
"""       
plot(   np.array([]),
        np.array([]),
        data_vol=V_TI1,
        data=(E_TI1-F_fit_TI1)-np.array(E_TI1-F_fit_TI1)[20],
        col=col6,
        pointtype="o",
        label="TI$\mathregular{^{full}}$")
"""        
p_temp = np.polyfit(V_fit_TI1,F_fit_TI1,7)   
F_temp = np.polyval(p_temp,V_TI1)       
TS_TI1 = (E_TI1-F_temp)-np.array(E_TI1-F_temp)[20]
plt.plot(   V_TI1/angstrom**3,
            (TS_TI1)/kjmol,
            color=col6,
            linestyle="",
            marker="o",
            label="TI$\mathregular{^{full}}$",
            markersize=4.)       

p_temp = np.polyfit(V_fit_TI6,F_fit_TI6,7)   
F_temp = np.polyval(p_temp,V_TI6)   
TS_TI6=(E_TI6-F_temp)-np.array(E_TI6-F_temp)[0]     
plt.plot(   V_TI6/angstrom**3,
            (TS_TI6)/kjmol,
            color=col5,
            linestyle="",
            marker="o",
            label="TI$\mathregular{^{ipol}}$",
            markersize=4)     
        
p_temp = np.polyfit(V_fit_TI7,F_fit_TI7,7)   
F_temp = np.polyval(p_temp,V_TI7)
TS_TI7=(E_TI7-F_temp)-np.array(E_TI7-F_temp)[0]        
plt.plot(   V_TI7/angstrom**3,
        (TS_TI7)/kjmol,
        color=col4,
        marker="o",
        linestyle="",
        label="TI$\mathregular{^{ipol2}}$",
        markersize=4)        
        
p_temp = np.polyfit(V_fit_TI4,F_fit_TI4,7)   
F_temp = np.polyval(p_temp,V_TI4)    
TS_TI4=(E_TI4-F_temp)-np.array(E_TI4-F_temp)[20]    
plt.plot(   V_TI4/angstrom**3,
        (TS_TI4)/kjmol,
        color=col3,
        marker="o",
        linestyle="",
        label="TI$\mathregular{^{mean}}$",
        markersize=4)         
        
"""
indexes = np.searchsorted(V_fit_TI3,V_TI3)
        
plot(   V_TI3[:-1],
        (E_TI3[:-1]-F_fit_TI3[indexes[:-1]])-np.min(E_TI3[:-1]-F_fit_TI3[indexes[:-1]]), 
        "m",
        linetype="--",
        label="TI$\mathregular{^{ipol}}$") 

indexes = np.searchsorted(V_fit_TI5,V_TI5)
        
plot(   V_TI5[:-1],
        (E_TI5[:-1]-F_fit_TI5[indexes[:-1]])-np.min(E_TI5[:-1]-F_fit_TI5[indexes[:-1]]), 
        "m",
        linetype="-.",
        label="TI$\mathregular{^{ipol2}}$") 
        
plot(   V_TI5[:-1],
        (E_TI5[:-1])-np.min(E_TI5[:-1]), 
        "g",
        linetype="-.",
        label="TI$\mathregular{^{ipol2}}$")                
"""

"""        
plot(   V_NMA1_242,
        TS_NMA1_242-np.array(TS_NMA1_242)[3],
        "r",
        linetype="--",
        label="2x4x2 NMA$\mathregular{^{mean}}$")
"""
plt.ylabel("$T \Delta S$ (kJ/mol)")    
plt.legend( loc="upper center",
            numpoints=1,
            ncol=2,bbox_to_anchor=(0.5, 1.3))
plt.ylim(-12,42)
plt.yticks([-10,0,10,20,30,40,50,60])
plt.savefig("S.pdf")
plt.savefig("S.png")
plt.clf()



####################################################################################

for i in xrange(3):
    for j in xrange(3):
        plotsettings()
        plt.plot(V_TI4/angstrom**3,RVECS_TI4[:,i,j]/angstrom,color=col3,marker=".",linestyle="")
        #plt.plot(V_TI2/angstrom**3,RVECS_TI2[:,i,j]/angstrom,color=col2,marker=".",linestyle="")
        plt.plot(V_TI6/angstrom**3,RVECS_TI6[:,i,j]/angstrom,color=col5,marker=".",linestyle="")
        plt.plot(V_TI7/angstrom**3,RVECS_TI7[:,i,j]/angstrom,color=col4,marker=".",linestyle="")
        #plt.plot(V_TI5_242/angstrom**3,RVECS_TI5_242[:,i,j]/angstrom,color=col5,marker="s",linestyle="")
        #plt.plot(V_TI4_242/angstrom**3,RVECS_TI4_242[:,i,j]/angstrom,color=col3,marker="s",linestyle="")
        plt.ylabel("$h_{" + str(i) + str(j) + "}$ ($\AA ^3$)")    
        plt.legend( loc="lower right",
                    handlelength=3,
                    ncol=2)
        plt.savefig("h/"+str(i)+str(j)+".png",bbox_inches="tight")
        plt.clf()

####################################################################################
titles=["a","b","c","alpha","beta","gamma"]
ylim=[[15,20],[12,17],[5,15]]
plotsettings()
for i in xrange(3):
        if i==0:
            plt.plot(V_TI2/angstrom**3,CELLS_TI2[:,i]/angstrom,color=col1,marker=".",linestyle="",markersize=4,label="TI$\mathregular{^{snap}}$")
            plt.plot(V_TI4/angstrom**3,CELLS_TI4[:,i]/angstrom,color=col3,marker=".",linestyle="",markersize=4,label="TI$\mathregular{^{mean}}$")
            plt.plot(V_TI6/angstrom**3,CELLS_TI6[:,i]/angstrom,color=col5,linestyle="-",linewidth=2,label="TI$\mathregular{^{ipol}}$")
            plt.plot(V_TI7/angstrom**3,CELLS_TI7[:,i]/angstrom,color=col4,linestyle="-",linewidth=2,label="TI$\mathregular{^{ipol2}}$")
        else:
            plt.plot(V_TI2/angstrom**3,CELLS_TI2[:,i]/angstrom,color=col1,marker=".",linestyle="",markersize=4)
            plt.plot(V_TI4/angstrom**3,CELLS_TI4[:,i]/angstrom,color=col3,marker=".",linestyle="",markersize=4)
            plt.plot(V_TI6/angstrom**3,CELLS_TI6[:,i]/angstrom,color=col5,linestyle="-",linewidth=2)
            plt.plot(V_TI7/angstrom**3,CELLS_TI7[:,i]/angstrom,color=col4,linestyle="-",linewidth=2)
        #plt.ylabel(titles[i] + " ($\AA$)")  
        #plt.ylim(ylim[i][0],ylim[i][1])  


plt.ylabel("Cell vector size ($\AA$)")
plt.ylim(5,20)
plt.legend( loc="lower right",
            handlelength=3,
            ncol=2)
plt.savefig("h/cellvecs.pdf",bbox_inches="tight")
plt.clf()
        
for i in xrange(3,6):
        plotsettings()
        plt.plot(V_TI4/angstrom**3,CELLS_TI4[:,i],color=col3,marker=".",linestyle="",markersize=8)
        plt.plot(V_TI6/angstrom**3,CELLS_TI6[:,i],color=col5,marker=".",linestyle="",markersize=8)
        plt.plot(V_TI7/angstrom**3,CELLS_TI7[:,i],color=col4,marker=".",linestyle="",markersize=8)
        plt.ylabel(titles[i] + " ($^{\circ}$)")    
        plt.legend( loc="lower right",
                    handlelength=3,
                    ncol=2)
        plt.ylim(80,100)
        plt.yticks([85,90,95])
        plt.savefig("h/"+titles[i]+".png",bbox_inches="tight")
        plt.clf()


       
####################################################################################

for i in xrange(3):
    for j in xrange(3):
        plotsettings()
        plt.plot(V_TI4_242/angstrom**3,RVECS_TI4_242[:,i,j]/angstrom,color=col3,marker="s",linestyle="")
        plt.plot(V_TI5_242/angstrom**3,RVECS_TI5_242[:,i,j]/angstrom,color=col5,marker="s",linestyle="")
        plt.ylabel("$h_{" + str(i) + str(j) + "}$ ($\AA ^3$)")    
        plt.legend( loc="lower right",
                    handlelength=3,
                    ncol=2)
        plt.savefig("h242/"+str(i)+str(j)+".png",bbox_inches="tight")
        plt.clf()


####################################################################################
        
for i in xrange(3):
    for j in xrange(3):
        plotsettings()
        plt.plot(V_TI1/angstrom**3,PTENS_TI1[:,i,j]/pascal/10**6,color=col6,marker=".",linestyle="")
        plt.plot(V_TI4/angstrom**3,PTENS_TI4[:,i,j]/pascal/10**6,color=col3,marker=".",linestyle="")
        #plt.plot(V_TI2/angstrom**3,PTENS_TI2[:,i,j]/pascal/10**6,color=col2,marker=".",linestyle="")
        plt.plot(V_TI6/angstrom**3,PTENS_TI6[:,i,j]/pascal/10**6,color=col5,marker=".",linestyle="")
        plt.plot(V_TI7/angstrom**3,PTENS_TI7[:,i,j]/pascal/10**6,color=col4,marker=".",linestyle="")
        #plt.plot(V_TI5_242/angstrom**3,PTENS_TI5_242[:,i,j]/pascal/10**6,color=col5,marker="s",linestyle="")
        #plt.plot(V_TI4_242/angstrom**3,PTENS_TI4_242[:,i,j]/pascal/10**6,color=col3,marker="s",linestyle="")
        plt.ylabel("$\sigma_{" + str(i) + str(j) + "}$ (MPa)")    
        plt.legend( loc="lower right",
                    handlelength=3,
                    ncol=2)
        plt.savefig("stress/sigma_"+str(i)+str(j)+".png",bbox_inches="tight")
        plt.clf()
        
####################################################################################
        
for i in xrange(3):
    for j in xrange(3):
        plotsettings()
        plt.plot(V_TI1_242/angstrom**3,PTENS_TI1_242[:,i,j]/pascal/10**6,color=col6,marker="s",linestyle="")
        plt.plot(V_TI4_242/angstrom**3,PTENS_TI4_242[:,i,j]/pascal/10**6,color=col3,marker="s",linestyle="")
        plt.plot(V_TI5_242/angstrom**3,PTENS_TI5_242[:,i,j]/pascal/10**6,color=col5,marker="s",linestyle="")
        plt.ylabel("$\sigma_{" + str(i) + str(j) + "}$ (MPa)")    
        plt.legend( loc="lower right",
                    handlelength=3,
                    ncol=2)
        plt.savefig("stress242/sigma_"+str(i)+str(j)+".png",bbox_inches="tight")
        plt.clf()

####################################################################################

print ""
print "#"*40
print ""
print "Metadynamics"
print "V:", np.round(crit_meta[0]/angstrom**3,3), " A**3"
print "+-sigma:" , np.round(crit_meta[1]/angstrom**3,3), " A**3"
print "F:", np.round(crit_meta[2]/kjmol,2), " kJ/mol"
print "+-sigma:", np.round(crit_meta[3]/kjmol,2), " kJ/mol"
print ""
'''
print "Metadynamics 2x4x2"
print "V:", np.round(crit_meta_242[0]/angstrom**3,3), " A**3"
print "+-sigma:" , np.round(crit_meta_242[1]/angstrom**3,1), " A**3"
print "F:", np.round(crit_meta_242[2]/kjmol,1), " kJ/mol"
print "+-sigma:", np.round(crit_meta_242[3]/kjmol,2), " kJ/mol"
print ""
'''
print "TI^full"
print "V:", np.round(crit_TI1[0]/angstrom**3,3), " A**3"
print "+-sigma:" , np.round(crit_TI1[1]*3/angstrom**3,3), " A**3"
print "F:", np.round(crit_TI1[2]/kjmol,2), " kJ/mol"
print "+-sigma:", np.round(crit_TI1[3]/kjmol,2), " kJ/mol"
print "V_np_sim: ", V_TI1[20]/angstrom**3, "E_np: ", E_TI1[20]/kjmol
print "V_lp_sim: ", V_TI1[-21]/angstrom**3, "E_np: ", E_TI1[-21]/kjmol
print "dE: ", E_TI1[-21]/kjmol-E_TI1[20]/kjmol, " dS" , TS_TI1[-21]/kjmol-TS_TI1[20]/kjmol
print ""
'''
print "TI^full 2x4x2"
print "V:", np.round(crit_TI1_242[0]/angstrom**3,3), " A**3"
print "+-sigma:" , np.round(crit_TI1_242[1]/angstrom**3,1), " A**3"
print "F:", np.round(crit_TI1_242[2]/kjmol,1), " kJ/mol"
print "+-sigma:", np.round(crit_TI1_242[3]/kjmol,2), " kJ/mol"
print "V_np_sim: ", V_TI1_242[5]/angstrom**3, "E_np: ", E_TI1_242[5]/kjmol
print "V_lp_sim: ", V_TI1_242[-5]/angstrom**3, "E_np: ", E_TI1_242[-5]/kjmol
print "dE: ", E_TI1_242[-5]/kjmol-E_TI1_242[5]/kjmol, " dS" , TS_TI1_242[-5]/kjmol-TS_TI1_242[5]/kjmol
print ""
'''
print "TI^snap"
print "V:", np.round(crit_TI2[0]/angstrom**3,3), " A**3"
print "+-sigma:" , np.round(crit_TI2[1]/angstrom**3,3), " A**3"
print "F:", np.round(crit_TI2[2]/kjmol,2), " kJ/mol"
print "+-sigma:", np.round(crit_TI2[3]/kjmol,2), " kJ/mol"
print ""
print "TI^ipol"
print "F:", np.round((F_fit_TI6[-1]-F_fit_TI6[0])/kjmol,3), " kJ/mol"
print "+-sigma:", np.round(crit_TI6[3]/kjmol,3), " kJ/mol"
print "V_np_sim: ", V_TI6[0]/angstrom**3, "E_np: ", E_TI6[0]/kjmol
print "V_lp_sim: ", V_TI6[-1]/angstrom**3, "E_np: ", E_TI6[-1]/kjmol
print "dE: ", E_TI6[-1]/kjmol-E_TI6[0]/kjmol, " dS" , TS_TI6[-1]/kjmol-TS_TI6[0]/kjmol
print ""
print "TI^ipol2"
print "F:", np.round((F_fit_TI7[-1]-F_fit_TI7[0])/kjmol,3), " kJ/mol"
print "+-sigma:", np.round(crit_TI7[3]/kjmol,3), " kJ/mol"
print "V_np_sim: ", V_TI7[0]/angstrom**3, "E_np: ", E_TI7[0]/kjmol
print "V_lp_sim: ", V_TI7[-1]/angstrom**3, "E_np: ", E_TI7[-1]/kjmol
print "dE: ", E_TI7[-1]/kjmol-E_TI7[0]/kjmol, " dS" , TS_TI7[-1]/kjmol-TS_TI7[0]/kjmol
print ""
print "TI^mean"
print "V:", np.round(crit_TI4[0]/angstrom**3,3), " A**3"
print "+-sigma:" , np.round(crit_TI4[1]/angstrom**3,3), " A**3"
print "F:", np.round(crit_TI4[2]/kjmol,2), " kJ/mol"
print "+-sigma:", np.round(crit_TI4[3]/kjmol,2), " kJ/mol"
print "V_np_sim: ", V_TI4[14]/angstrom**3, "E_np: ", E_TI4[14]/kjmol
print "V_lp_sim: ", V_TI4[-20]/angstrom**3, "E_np: ", E_TI4[-20]/kjmol
print "dE: ", E_TI4[-20]/kjmol-E_TI4[14]/kjmol, " dS:" , TS_TI4[-20]/kjmol-TS_TI4[14]/kjmol
print ""
'''
print "TI^mean 2x4x2"
print "V:", np.round(crit_TI4_242[0]/angstrom**3,3), " A**3"
print "+-sigma:" , np.round(crit_TI4_242[1]/angstrom**3,1), " A**3"
print "F:", np.round(crit_TI4_242[2]/kjmol,1), " kJ/mol"
print "+-sigma:", np.round(crit_TI4_242[3]/kjmol,2), " kJ/mol"
print "V_np_sim: ", V_TI4_242[5]/angstrom**3, "E_np: ", E_TI4_242[5]/kjmol
print "V_lp_sim: ", V_TI4_242[-5]/angstrom**3, "E_np: ", E_TI4_242[-5]/kjmol
print "dE: ", E_TI4_242[-5]/kjmol-E_TI4_242[5]/kjmol, " dS:" , TS_TI4_242[-5]/kjmol-TS_TI4_242[5]/kjmol
print ""
print "TI^ipol 2x4x2"
print "F:", np.round((F_fit_TI5_242[-1]-F_fit_TI5_242[0])/kjmol,1), " kJ/mol"
print "+-sigma:", np.round(crit_TI5_242[3]/kjmol,2), " kJ/mol"
print "V_np_sim: ", V_TI5_242[0]/angstrom**3, "E_np: ", E_TI5_242[0]/kjmol
print "V_lp_sim: ", V_TI5_242[-1]/angstrom**3, "E_np: ", E_TI5_242[-1]/kjmol
print "dE: ", E_TI5_242[-1]/kjmol-E_TI5_242[0]/kjmol, " dS:" , TS_TI5_242[-1]/kjmol-TS_TI5_242[0]/kjmol
print ""
'''
####################################################################################


print "#"*40
print " DUT49"
print "#"*40
print ""
print "MTD NVsT"
print "V:", np.round(np.real(crit_MTD1_DUT[0])/angstrom**3,0), " A**3"
print "+-sigma:" , np.round(np.real(crit_MTD1_DUT[1])/angstrom**3,0), " A**3"
print ""
print "MTD NVhT"
print "V:", np.round(np.real(crit_MTD2_DUT[0])/angstrom**3,0), " A**3"
print "+-sigma:" , np.round(np.real(crit_MTD2_DUT[1])/angstrom**3,0), " A**3"
print ""
print "TI^full_2 (OP)"
print "V:", np.round(np.real(crit_TI1_DUT[0])/angstrom**3,0), " A**3"
print "+-sigma:" , np.round(np.real(crit_TI1_DUT[1])/angstrom**3,0), " A**3"
print "V_E:", np.round(V_TI1_DUT[93]/angstrom**3)
print "E: " , np.round(E_TI1_DUT[93]/kjmol)
print ""
print "TI^resc_2 (OP)"
print "V:", np.round(np.real(crit_TI2_DUT[0])/angstrom**3,0), " A**3"
print "+-sigma:" , np.round(np.real(crit_TI2_DUT[1])/angstrom**3,0), " A**3"
print "V_E:", np.round(V_TI2_DUT[93]/angstrom**3)
print "E: " , np.round(E_TI2_DUT[93]/kjmol)
print ""
print "TI^full_1 (CP)"
print "F:", np.round(crit_TI3_DUT[2]/kjmol,1), " kJ/mol"
print "+-sigma:", np.round(crit_TI3_DUT[3]/kjmol,2), " kJ/mol"
print "V:", np.round(crit_TI3_DUT[0]/angstrom**3,0), " A**3"
print "+-sigma:" , np.round(crit_TI3_DUT[1]/angstrom**3,0), " A**3"
print "V_E:", np.round(V_TI3_DUT[93]/angstrom**3),np.round(V_TI3_DUT[21]/angstrom**3),np.round(V_TI3_DUT[10]/angstrom**3)
print "E:",  np.round((E_TI3_DUT[93])/kjmol) ,  np.round((E_TI3_DUT[10])/kjmol) , "dS: ", np.round((E_TI3_DUT[93]-E_TI3_DUT[10] - crit_TI3_DUT[2,1])/kjmol)
print ""
print "TI^resc_1 (CP)"
print "F:", np.round(crit_TI4_DUT[2]/kjmol,1), " kJ/mol"
print "+-sigma:", np.round(crit_TI4_DUT[3]/kjmol,2), " kJ/mol"
print "V:", np.round(crit_TI4_DUT[0]/angstrom**3,0), " A**3"
print "+-sigma:" , np.round(crit_TI4_DUT[1]/angstrom**3,0), " A**3"
print "V_E:", np.round(V_TI4_DUT[93]/angstrom**3),np.round(V_TI4_DUT[21]/angstrom**3),np.round(V_TI4_DUT[10]/angstrom**3)
print "dE:",  np.round((E_TI4_DUT[93])/kjmol) ,  np.round((E_TI4_DUT[10])/kjmol) , "dS: ", np.round((E_TI4_DUT[93]-E_TI4_DUT[10] - crit_TI4_DUT[2,1])/kjmol)
print ""
print "#"*40
print ""
####################################################################################

#                 DUT49


f,ax=plt.subplots()

plotsettings()
plot_fill_P(ax,V_fit_TI1_DUT,P_fit_TI1_DUT,label="TI$\mathregular{^{full}}\mathregular{_{2}}$",lt="-",col=col1,alpha=0.5)
plot_fill_P(ax,V_fit_TI2_DUT,P_fit_TI2_DUT,label="TI$\mathregular{^{resc}}\mathregular{_{2}}$",lt="-",col=col2,alpha=0.5)
plot_fill_P(ax,V_fit_TI3_DUT,P_fit_TI3_DUT,label="TI$\mathregular{^{full}}\mathregular{_{1}}$",lt="-",col=col3,alpha=0.5)
plot_fill_P(ax,V_fit_TI4_DUT,P_fit_TI4_DUT,label="TI$\mathregular{^{resc}}\mathregular{_{1}}$",lt="-",col=col5,alpha=0.5)


plt.plot(V_TI1_DUT/angstrom**3,P_TI1_DUT/pascal/10**6,"o",linestyle="",markerfacecolor=col1,markeredgecolor="k",markersize=3.,alpha=0.5)

plt.plot(V_TI2_DUT/angstrom**3,P_TI2_DUT/pascal/10**6,"o",linestyle="",markerfacecolor=col2,markeredgecolor="k",markersize=3.,alpha=0.5)

plt.plot(V_TI3_DUT/angstrom**3,P_TI3_DUT/pascal/10**6,"o",linestyle="",markerfacecolor=col3,markeredgecolor="k",markersize=3.,alpha=0.5)

plt.plot(V_TI4_DUT/angstrom**3,P_TI4_DUT/pascal/10**6,"o",linestyle="",markerfacecolor=col5,markeredgecolor="k",markersize=3.,alpha=0.5)


#plot_fill_P(ax,V_fit_TI1_242,P_fit_TI1_242,label="$\mathregular{TI^{full}_{2x4x2}}$",lt="--",col="b",alpha=0.25)

#plot_fill_P(ax,V_fit_TI4_242,P_fit_TI4_242,label="$\mathregular{TI^{mean}_{2x4x2}}$",lt="--",col="m",alpha=0.4)

plt.ylabel("$P$ (MPa)")
plt.ylim(-0.2e3,0.2e3)
plt.yticks([-0.1e3,0,0.1e3]) 
#plt.title("1x2x1 vs 2x4x2 unit cell")
plt.xlabel("$V$ ($\AA ^3$)")
plt.xlim(40000,100000)
plt.xticks([40000,60000,80000,100000])
plt.legend(loc="upper center",ncol=2,numpoints=1,bbox_to_anchor=(0.5, 1.3))
#plt.tight_layout()
plt.savefig("P_DUT.png")
plt.savefig("P_DUT.pdf")
#plt.show()
plt.clf()


####################################################################################

plotsettings()       

plot(   V_fit_TI1_DUT,
        F_fit_TI1_DUT  ,
        col1,
        crit_points=crit_TI1_DUT,
        crit_points_npy=-100,
        linetype="-",
        #data_vol=V_meta,
        #data=F_meta,
        every=5,
        label="TI$\mathregular{^{full}}\mathregular{_{2}}$")
        
plot(   V_fit_TI2_DUT,
        F_fit_TI2_DUT  ,
        col2,
        crit_points=crit_TI2_DUT,
        crit_points_npy=-100,
        linetype="-",
        #data_vol=V_meta,
        #data=F_meta,
        every=5,
        label="TI$\mathregular{^{resc}}\mathregular{_{2}}$")
        
plot(   V_fit_TI3_DUT,
        F_fit_TI3_DUT  ,
        col3,
        crit_points=crit_TI3_DUT,
        crit_points_npy=-100,
        linetype="-",
        #data_vol=V_meta,
        #data=F_meta,
        every=5,
        label="TI$\mathregular{^{full}}\mathregular{_{1}}$")
        
plot(   V_fit_TI4_DUT,
        F_fit_TI4_DUT,
        col5,
        linetype="-",
        #crit_points=crit_TI3,
        crit_points_npy=-100,
        #dashes=(10,3,3,3),
        label="TI$\mathregular{^{resc}}\mathregular{_{1}}$")
     

      
plt.ylabel("$F$ (kJ/mol)")
#plt.ylim(-10,50)
#plt.yticks([0,10,20,30,40]) 
#plt.title("TI$\mathregular{^{full}}$ gives the only precise estimate of the free energy")        
plt.legend( loc="upper center",
            handlelength=3,
            ncol=2,bbox_to_anchor=(0.5, 1.3))
#plt.tight_layout()
plt.xlim(40000,100000)
plt.xticks([40000,60000,80000,100000])
plt.ylim(0,2000)
plt.savefig("F_DUT49.png")
plt.savefig("F_DUT49.pdf")
#plt.show()    
plt.clf() 

####################################################################################

plotsettings()       
        
plot(   V_fit_TI3_DUT,
        F_fit_TI3_DUT  ,
        col3,
        crit_points=crit_TI3_DUT,
        crit_points_npy=-100,
        linetype="-",
        #data_vol=V_meta,
        #data=F_meta,
        every=5,
        label="TI$\mathregular{^{full}}\mathregular{_{1}}$")
        
plot(   V_fit_TI4_DUT,
        F_fit_TI4_DUT,
        col5,
        linetype="-",
        #crit_points=crit_TI3,
        crit_points_npy=-100,
        #dashes=(10,3,3,3),
        label="TI$\mathregular{^{resc}}\mathregular{_{1}}$")
     

      
plt.ylabel("$F$ (kJ/mol)")
#plt.ylim(-10,50)
#plt.yticks([0,10,20,30,40]) 
#plt.title("TI$\mathregular{^{full}}$ gives the only precise estimate of the free energy")        
plt.legend( loc="upper center",
            handlelength=3,
            ncol=2,bbox_to_anchor=(0.5, 1.3))
#plt.tight_layout()
plt.xlim(44600,47000)
plt.xticks([44600,45200,45800,46400,47000])
plt.ylim(1055,1095)
plt.yticks([1060,1070,1080,1090])
plt.savefig("F_DUT49_zoomed.png")
plt.savefig("F_DUT49_zoomed.pdf")
#plt.show()    
plt.clf() 

####################################################################################

plotsettings()       

plt.plot(   V_TI1_DUT/angstrom**3,
        E_TI1_DUT/kjmol  ,
        color=col1,
        linestyle="",
        marker="o",
        markersize=3.,
        label="TI$\mathregular{^{full}}\mathregular{_{2}}$")
        
plt.plot(   V_TI2_DUT/angstrom**3,
        E_TI2_DUT/kjmol  ,
        color=col2,
        linestyle="",
        marker="o",
        markersize=3.,
        label="TI$\mathregular{^{resc}}\mathregular{_{2}}$")
        
plt.plot(   V_TI3_DUT/angstrom**3,
        E_TI3_DUT/kjmol  ,
        color=col3,
        linestyle="",
        marker="o",
        markersize=3.,
        label="TI$\mathregular{^{full}}\mathregular{_{1}}$")
        
plt.plot(   V_TI4_DUT/angstrom**3,
        E_TI4_DUT/kjmol,
        color=col5,
        linestyle="",
        marker="o",
        markersize=3.,
        label="TI$\mathregular{^{resc}}\mathregular{_{1}}$")
     

      
plt.ylabel("$E$ (kJ/mol)")
plt.legend( loc="upper center",
            handlelength=3,
            ncol=2,bbox_to_anchor=(0.5, 1.3),numpoints=1
            )
plt.xlim(40000,100000)
plt.xticks([40000,60000,80000,100000])
plt.savefig("E_DUT49.png")
plt.savefig("E_DUT49.pdf")
#plt.show()    
plt.clf() 





























'''
####################################################################################

assert False

#FML 7/4/17
F
plotsettings()
        
        
plt.plot(V_TI1_242/angstrom**3,P_TI1_242/pascal/10**6,"bx",label="TI method 1",markersize=15,markeredgewidth=2)
plt.plot(V_TI4_242/angstrom**3,P_TI4_242/pascal/10**6,"rx",label="TI method 2",markersize=15,markeredgewidth=2)
             

plt.ylabel("$P$ (MPa)",fontsize=20)
plt.ylim(-0.75e3,0.75e3)

plt.xticks([750,1000,1250,1500],fontsize=20) 
plt.xlabel("$V$ ($\AA ^3$)",fontsize=20)
plt.yticks([-0.5e3,0,0.5e3],fontsize=20) 
#plt.title("1x2x1 vs 2x4x2 unit cell") 
plt.grid()
plt.legend(loc="upper center",ncol=2)
plt.savefig("FML/P_242_1.png")
plt.clf()    

plotsettings()
        
plt.plot(V_TI1_242/angstrom**3,P_TI1_242/pascal/10**6,"bx",markersize=15,markeredgewidth=2)
plt.plot(V_TI4_242/angstrom**3,P_TI4_242/pascal/10**6,"rx",markersize=15,markeredgewidth=2)
plt.plot(V_fit_TI1_242/angstrom**3,P_fit_TI1_242/pascal/10**6,"b-",label="TI method 1",lw=3)
plt.plot(V_fit_TI4_242/angstrom**3,P_fit_TI4_242/pascal/10**6,"r-",label="TI method 2",lw=3)
             

plt.ylabel("$P$ (MPa)",fontsize=20)
plt.ylim(-0.75e3,0.75e3)

#plt.xticks([750,1000,1250,1500],fontsize=20) 
#plt.xlabel("$V$ ($\AA ^3$)",fontsize=20)
plt.yticks([-0.5e3,0,0.5e3],fontsize=20) 
#plt.title("1x2x1 vs 2x4x2 unit cell") 
plt.grid()
plt.legend(loc="upper center",ncol=2)
plotsettings()
plt.savefig("FML/P_242_2.png")
#plt.show()   
plt.clf() 

plotsettings()

import matplotlib as mpl
mpl.rcParams['lines.linewidth'] = 2


plot(   V_fit_TI1_242,
        F_fit_TI1_242, 
        "b",
        #crit_points=crit_TI1_242,
        #crit_points_npy=-100,
        linetype="--",
        #data_vol=V_meta_242,
        #data=F_meta_242,
        every=5,
        label="TI method 1")
        
plot(   V_fit_TI4_242,
        F_fit_TI4_242, 
        "r",
        #crit_points=crit_TI4_242,
        #crit_points_npy=-100,
        linetype="--",
        #data_vol=V_meta_242,
        #data=F_meta_242,
        every=5,
        label="TI method 2")


plt.ylabel("$F$ (kJ/mol)",fontsize=20)
plt.xlabel("$V$ ($\AA ^3$)",fontsize=20)
plt.ylim(-2,42)
plt.yticks([0,20,40],fontsize=20) 

#plt.title("1x2x1 vs 2x4x2 unit cell") 
plt.grid()
plt.xticks([750,1000,1250,1500],fontsize=20) 
plt.legend(loc="upper center",ncol=2)
plt.savefig("FML/F_242_1.png")
plt.clf()   

plotsettings()

import matplotlib as mpl



plot(   V_fit_TI1_242,
        F_fit_TI1_242, 
        "b",
        crit_points=crit_TI1_242,
        crit_points_npy=-100,
        linetype="--",
        #data_vol=V_meta_242,
        #data=F_meta_242,
        every=5,
        label="TI method 1")
        
plot(   V_fit_TI4_242,
        F_fit_TI4_242, 
        "r",
        crit_points=crit_TI4_242,
        crit_points_npy=-100,
        linetype="--",
        #data_vol=V_meta_242,
        #data=F_meta_242,
        every=5,
        label="TI method 2")


plt.ylabel("$F$ (kJ/mol)",fontsize=20)
plt.xlabel("$V$ ($\AA ^3$)",fontsize=20)
plt.ylim(-2,42)
plt.yticks([0,20,40],fontsize=20) 

#plt.title("1x2x1 vs 2x4x2 unit cell") 
plt.grid()
plt.xticks([750,1000,1250,1500],fontsize=20) 
plt.legend(loc="upper center",ncol=2)
plt.savefig("FML/F_242_2.png")
plt.clf()    





#plt.ylim(-5,45)

'''
###################### Added Sven #######################

plotsettings()

f, (ax, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw = {'height_ratios':[1, 2]})

ax.plot(	V_fit_meta/angstrom**3,
			F_fit_meta/kjmol,
			color = col_MTD,
			linestyle = '-',
			linewidth = 3,
			label = "$\mathregular{MTD}$" )

ax.plot(	V_fit_TI1/angstrom**3,
			F_fit_TI1/kjmol,
			color = col_TIfull,
			linestyle = '-',
			label = "$\mathregular{TI_{full}}$")

ax.plot(	V_fit_TI3_old/angstrom**3,
			F_fit_TI3_old/kjmol,
			color = col_ipol1,
			linestyle = '-',
			label = "$\mathregular{TI^{ipol1}_{\mathbf{h}_0}}$")
			


ax2.plot(	V_fit_meta/angstrom**3,
			F_fit_meta/kjmol,
			color = col_MTD,
			linestyle = '-',
			linewidth = 3,
			label = "$\mathregular{MTD}$" )

ax2.plot(	V_fit_TI1/angstrom**3,
			F_fit_TI1/kjmol,
			color = col_TIfull,
			linestyle = '-',
			label = "$\mathregular{TI_{full}}$")

ax2.plot(	V_fit_TI3_old/angstrom**3,
			F_fit_TI3_old/kjmol,
			color = col_ipol1,
			linestyle = '-',
			label = "$\mathregular{TI^{ipol1}_{\mathbf{h}_0}}$")

ax.set_ylim(180, 200)
ax2.set_ylim(0,40)

ax.set_yticks([180, 190, 200])
ax2.set_yticks([0, 10, 20, 30, 40])

ax.spines['bottom'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax.xaxis.tick_top()
ax.tick_params(labeltop='off')
ax2.xaxis.tick_bottom()

d = .015
kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
ax.plot((-d, +d), (-d, +d), **kwargs)
ax.plot((1-d, 1+d), (-d, +d), **kwargs)

kwargs.update(transform=ax2.transAxes)
ax2.plot((-d, +d), (1-d, 1+d), **kwargs)
ax2.plot((1-d, 1+d), (1-d, 1+d), **kwargs)


'''
plot(   V_fit_meta,
        F_fit_meta,
        col_MTD,
        crit_points=crit_meta,
        crit_points_npy=-100,
        linetype="-",
        every=5,
        label="$\mathregular{MTD}$")
        
plot(   V_fit_TI1,
        F_fit_TI1,
		col_TIfull,
        crit_points=crit_TI1,
        crit_points_npy=-100,
        linetype="-",
        label="$\mathregular{TI_{full}}$")

plot(   V_fit_TI6,
        F_fit_TI6,
		col_ipol1,
        crit_points=crit_TI6,
        crit_points_npy=-100,
        linetype="-",
        label="$\mathregular{TI^{ipol1}_{\mathbf{h}_0}}$")


plot(   V_fit_TI3_old,
        F_fit_TI3_old,
		col_ipol1,
        crit_points=crit_TI3_old,
        crit_points_npy=-100,
        linetype="-",
        label="$\mathregular{TI^{ipol1}_{\mathbf{h}_0}}$")
'''

plt.xlabel("Unit cell volume $V$ [$\mathrm{\AA}^3$]")
plt.ylabel("Free energy $F$ [kJ/mol]")
#plt.yticks([0,10,20,30,40]) 
    
plt.legend( loc="upper center",
            handlelength=2,
            ncol=3,
            borderpad=0.2,
            borderaxespad=0.25)

plt.savefig("Figure3a.png", bbox_inches = 'tight')
plt.savefig("Figure3a.pdf", format='pdf', bbox_inches = 'tight')
plt.savefig("Figure3a.svg", bbox_inches = 'tight')

plt.clf()
      
plot(   V_fit_TI3_DUT,
        F_fit_TI3_DUT,
		col_TIfull,
        crit_points=crit_TI3_DUT,
        crit_points_npy=-100,
        linetype="-",
        every=5,
        label="$\mathregular{TI_{full}}$")
        
plot(   V_fit_TI4_DUT,
        F_fit_TI4_DUT,
		col_ipol1,
        linetype="-",
        crit_points_npy=-100,
        label="$\mathregular{TI^{ipol}_{\mathbf{h}_0}}$")
     

plt.xlabel("Unit cell volume $V$ [$\mathrm{\AA}^3$]")
plt.ylabel("Free energy $F$ [kJ/mol]")      
plt.legend( loc="upper center",
            handlelength=3,
            ncol=2)

plt.xlim(40000,100000)
plt.xticks([40000,50000,60000,70000,80000,90000,100000])
plt.ylim(0,1500)


plt.savefig("Figure3b.png", bbox_inches = 'tight')
plt.savefig("Figure3b.pdf", format='pdf', bbox_inches = 'tight')
plt.savefig("Figure3b.svg", bbox_inches = 'tight')


plt.clf()


#############

plt.subplot(211)

# TI_full
plt.plot(	V_TI1/angstrom**3,
			P_TI1/pascal/10**6,
			"o",
			linestyle="",
			markerfacecolor=col_TIfull,
			markeredgecolor="k",
			markersize=4,
			label = "$\mathregular{TI_{full}}$")

# TI_snap
plt.plot(	V_TI2/angstrom**3,
			P_TI2/pascal/10**6,
			"o",
			linestyle="",
			markerfacecolor=col_snap,
			markeredgecolor="k",
			markersize=4,
			label = "$\mathregular{TI^{snap}_{\mathbf{h}_0}}$")

# TI_mean
plt.plot(	V_TI4/angstrom**3,
			P_TI4/pascal/10**6,
			"o",
			linestyle="",
			markerfacecolor=col_mean,
			markeredgecolor="k",
			markersize=4,
			label = "$\mathregular{TI^{mean}_{\mathbf{h}_0}}$")

# TI_ipol1
plt.plot(	V_TI6/angstrom**3,
			P_TI6/pascal/10**6,
			"o",
			linestyle="",
			markerfacecolor=col_ipol1,
			markeredgecolor="k",
			markersize=4,
			label = "$\mathregular{TI^{ipol1}_{\mathbf{h}_0}}$")


# TI_ipol2
plt.plot(	V_TI7/angstrom**3,
			P_TI7/pascal/10**6,
			"o",
			linestyle="",
			markerfacecolor=col_ipol2,
			markeredgecolor="k",
			markersize=4,
			label = "$\mathregular{TI^{ipol2}_{\mathbf{h}_0}}$")

plt.ylim(-1000, 1500)
plt.xlabel("Unit cell volume $V$ [$\mathrm{\AA}^3$]")
plt.ylabel("Hydrostatic pressure $P$ [MPa]")
    
plt.legend( loc="upper center",
            handlelength=2,
            ncol=3,
            borderpad=0.2,
            borderaxespad=0.25)

plt.subplot(212)

# TI_full
plt.plot(	V_TI1/angstrom**3,
			P_TI1/pascal/10**6,
			"s",
			linestyle="",
			markerfacecolor=col_TIfull,
			markeredgecolor="k",
			markersize = 3,
			label = "$\mathregular{TI_{full}}$")

# TI_ipol1

plt.plot(	V_TI6/angstrom**3,
			P_ATI_TI6/pascal/10**6,
			"s",
			linestyle="",
			markerfacecolor=col_ipol1,
			markeredgecolor="k",
			markersize = 3,
			label = "$\mathregular{TI^{ipol1}_{\mathbf{h}_0}}$")

# TI_ipol2

plt.plot(	V_TI7/angstrom**3,
			P_ATI_TI7/pascal/10**6,
			"s",
			linestyle="",
			markerfacecolor=col_ipol2,
			markeredgecolor="k",
			markersize = 3,
			label = "$\mathregular{TI^{ipol2}_{\mathbf{h}_0}}$")

plt.xlabel("Unit cell volume $V$ [$\mathrm{\AA}^3$]")
plt.ylabel("Generalized pressure $P_a$ [MPa]")
    
plt.legend( loc="upper center",
            handlelength=2,
            ncol=3,
            borderpad=0.2,
            borderaxespad=0.25)

fig = plt.gcf()
fig.set_size_inches(6.5, 6)
plt.tight_layout()
plt.savefig("Figure4.png", bbox_inches = 'tight')
plt.savefig("Figure4.pdf", format='pdf', bbox_inches = 'tight')
plt.savefig("Figure4.svg", bbox_inches = 'tight')


######

plotsettings()
        

#TI_full

plt.plot(   V_fit_TI1/angstrom**3,
        	F_fit_TI1/kjmol,
        	color = col_TIfull,
			linestyle = '-',
			label = "$\mathregular{TI_{full}}$")

# TI_ipol1
plt.plot(	V_fit_TI6/angstrom**3,
			F_fit_TI6/kjmol,
			color=col_ipol1,
			linestyle = '-',
			label = "$\mathregular{TI^{ipol1}_{\mathbf{h}_0}}$")


# TI_ipol2
plt.plot(	V_fit_TI7/angstrom**3,
			F_fit_TI7/kjmol,
			color=col_ipol2,
			linestyle = '-',
			label = "$\mathregular{TI^{ipol2}_{\mathbf{h}_0}}$")
      
plt.xlabel("Unit cell volume $V$ [$\mathrm{\AA}^3$]")
plt.ylabel("Free energy $F$ [kJ/mol]")
plt.ylim(0,40)
plt.xlim(700,1600)
plt.yticks([0,10,20,30,40])
plt.xticks([700,800,900,1000,1100,1200,1300,1400,1500])
plt.legend( loc="upper center",
            handlelength=3,
            ncol=3,bbox_to_anchor=(0.5, 0.9))

plt.savefig("Figure5.png", bbox_inches = 'tight')
plt.savefig("Figure5.pdf", format='pdf', bbox_inches = 'tight')
plt.savefig("Figure5.svg", bbox_inches = 'tight')

plt.clf() 


###


plotsettings()


plt.subplot(211)

E0 = np.amin(E_TI1)

print('MIL-53')
print(E0/kjmol)

## TI_full
plt.plot(   V_TI1/angstrom**3,
        	(E_TI1-E0)/kjmol,
        	color = col_TIfull,
			linestyle = '-',
			label = "$\mathregular{TI_{full}}$")
        
        
# TI_ipol1
plt.plot(	V_TI6/angstrom**3,
			(E_TI6-E0)/kjmol,
			color=col_ipol1,
			linestyle = '-',
			label = "$\mathregular{TI^{ipol1}_{\mathbf{h}_0}}$")


# TI_ipol2
plt.plot(	V_TI7/angstrom**3,
			(E_TI7-E0)/kjmol,
			color=col_ipol2,
			linestyle = '-',
			label = "$\mathregular{TI^{ipol2}_{\mathbf{h}_0}}$")   

plt.xlabel('Unit cell volume $V$ [$\mathrm{\AA}^3$]')
plt.ylabel('Energy $\Delta E$ [kJ/mol]')

plt.legend( loc="upper center",
            numpoints=1,
            ncol=2,
            bbox_to_anchor=(0.5, 0.95))

plt.subplot(212)

## TI_full

p_temp = np.polyfit(V_fit_TI1,F_fit_TI1,7)   
F_temp = np.polyval(p_temp,V_TI1)       
TS_TI1 = (E_TI1-F_temp)-np.array(E_TI1-F_temp)[20]

plt.plot(   V_TI1/angstrom**3,
            (TS_TI1)/kjmol,
            color = col_TIfull,
			linestyle = '-',
			label = "$\mathregular{TI_{full}}$")

## TI_ipol1   

p_temp = np.polyfit(V_fit_TI6,F_fit_TI6,7)   
F_temp = np.polyval(p_temp,V_TI6)   
TS_TI6=(E_TI6-F_temp)-np.array(E_TI6-F_temp)[0]

plt.plot(   V_TI6/angstrom**3,
            (TS_TI6)/kjmol,
            color=col_ipol1,
			linestyle = '-',
			label = "$\mathregular{TI^{ipol1}_{\mathbf{h}_0}}$")

## TI_ipol2
        
p_temp = np.polyfit(V_fit_TI7,F_fit_TI7,7)   
F_temp = np.polyval(p_temp,V_TI7)
TS_TI7=(E_TI7-F_temp)-np.array(E_TI7-F_temp)[0]
     
plt.plot(   V_TI7/angstrom**3,
        	(TS_TI7)/kjmol,
        	color=col_ipol2,
			linestyle = '-',
			label = "$\mathregular{TI^{ipol2}_{\mathbf{h}_0}}$")    
        
plt.xlabel('Unit cell volume $V$ [$\mathrm{\AA}^3$]')
plt.ylabel('Entropy $T \Delta S$ [kJ/mol]')

plt.legend( loc="upper center",
            numpoints=1,
            ncol=2,
            bbox_to_anchor=(0.5, 0.95))

fig = plt.gcf()
fig.set_size_inches(4, 7)
plt.tight_layout()
plt.savefig("Figure6.png", bbox_inches = 'tight')
plt.savefig("Figure6.pdf", format='pdf', bbox_inches = 'tight')
plt.savefig("Figure6.svg", bbox_inches = 'tight')

plt.clf()




########

plotsettings()


plt.subplot(211)

E0_DUT = np.amin(E_TI3_DUT)


print('DUT-49')
print(E0_DUT/kjmol)

## TI_full
plt.plot(   V_TI3_DUT/angstrom**3,
        	(E_TI3_DUT-E0_DUT)/kjmol  ,
        	color = col_TIfull,
			linestyle = '-',
			label = "$\mathregular{TI_{full}}$")

# TI_ipol        
plt.plot(   V_TI4_DUT/angstrom**3,
        	(E_TI4_DUT-E0_DUT)/kjmol,
        	color=col_ipol1,
			linestyle = '-',
			label = "$\mathregular{TI^{ipol}_{\mathbf{h}_0}}$")


plt.xlabel('Unit cell volume $V$ [$\mathrm{\AA}^3$]')
plt.ylabel('Energy $\Delta E$ [kJ/mol]')

plt.legend( loc="upper center",
            numpoints=1,
            ncol=2,
            bbox_to_anchor=(0.5, 0.95))

plt.subplot(212)

## TI_full

p_temp = np.polyfit(V_fit_TI3_DUT,F_fit_TI3_DUT,7)   
F_temp = np.polyval(p_temp,V_TI3_DUT)       

abs_min = min(P_fit_TI3_DUT[:len(P_fit_TI3_DUT)/5], key=abs)
ind_fit_min = np.where(P_fit_TI3_DUT[:len(P_fit_TI3_DUT)/5] == abs_min)[0]
vol_CP = V_fit_TI3_DUT[ind_fit_min]
ind_min = np.argmin((V_TI3_DUT - vol_CP)**2)

TS_TI3_DUT = (E_TI3_DUT-F_temp)-np.array(E_TI3_DUT-F_temp)[ind_min]

plt.plot(   V_TI3_DUT/angstrom**3,
            (TS_TI3_DUT)/kjmol,
            color = col_TIfull,
			linestyle = '-',
			label = "$\mathregular{TI_{full}}$")

## TI_ipol

p_temp = np.polyfit(V_fit_TI4_DUT,F_fit_TI4_DUT,7)   
F_temp = np.polyval(p_temp,V_TI4_DUT)       

abs_min = min(P_fit_TI4_DUT[:len(P_fit_TI4_DUT)/7], key=abs)
ind_fit_min = np.where(P_fit_TI4_DUT[:len(P_fit_TI4_DUT)/7] == abs_min)[0]
vol_CP = V_fit_TI4_DUT[ind_fit_min]
ind_min = np.argmin((V_TI4_DUT - vol_CP)**2)

TS_TI4_DUT = (E_TI4_DUT-F_temp)-np.array(E_TI4_DUT-F_temp)[ind_min]

plt.plot(   V_TI4_DUT/angstrom**3,
            (TS_TI4_DUT)/kjmol,
            color=col_ipol1,
			linestyle = '-',
			label = "$\mathregular{TI^{ipol}_{\mathbf{h}_0}}$")
        
plt.xlabel('Unit cell volume $V$ [$\mathrm{\AA}^3$]')
plt.ylabel('Entropy $T \Delta S$ [kJ/mol]')

plt.legend( loc="upper center",
            numpoints=1,
            ncol=2,
            bbox_to_anchor=(0.5, 0.95))

fig = plt.gcf()
fig.set_size_inches(4, 7)
plt.tight_layout()
plt.savefig("Figure7.png", bbox_inches = 'tight')
plt.savefig("Figure7.pdf", format='pdf', bbox_inches = 'tight')
plt.savefig("Figure7.svg", bbox_inches = 'tight')

plt.clf()

#########



plotsettings()


plt.plot(	V_TI1/angstrom**3,
			P_TI1/(1e6*pascal),
			color = col_TIfull,
			linestyle = '-',
			label = r'1 $\times$ 2 $\times$ 1')

plt.plot(	V_TI1_242/angstrom**3,
			P_TI1_242/(1e6*pascal),
			color = col_MTD,
			linestyle = '-',
			label = r'2 $\times$ 4 $\times$ 2')
        
plt.xlabel('Unit cell volume $V$ [$\mathrm{\AA}^3$]')
plt.ylabel('Pressure $P$ [MPa]')

plt.legend( loc="upper center",
            numpoints=1,
            ncol=2,
            bbox_to_anchor=(0.5, 0.95))

fig = plt.gcf()
plt.savefig("FigureS5.png", bbox_inches = 'tight')
plt.savefig("FigureS5.pdf", format='pdf', bbox_inches = 'tight')
plt.savefig("FigureS5.svg", bbox_inches = 'tight')

plt.clf()
