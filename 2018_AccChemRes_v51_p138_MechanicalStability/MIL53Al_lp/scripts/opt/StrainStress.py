''' Stel je vertrekt van een sample (strain en stress vector gekend) ==> Hoe wordt matrix dan gevormd? '''
''' Cij vector gedefinieerd als kolomvector met elementen:
 (C11 C21 C31 C41 C51 C61 C22 C32 C42 C52 C62 C33 C43 C53 C63 C44 C54 C64 C55 C65 C66)^T
 Voigt notation
 Met onderstaande keuze van matrix
  C' = cell tensor = C*epsilon 
   = ( a1 a2 a3 ) * ( 1+e1 e6/2 e5/2 )
     ( b1 b2 b3 ) * ( e6/2 1+e2 e4/2 )
     ( c1 c2 c3 ) * ( e5/2 e4/2 1+e3 )
  
  strain = epsilon = (e1 e2 e3 e4 e5 e6)^T
  stress = sigma = (s1 s2 s3 s4 s5 s6)^T
  
  De stress wordt voor verschillende strains berekend en vervolgens een least-square fit om de elastische constanten te bepalen.
  
'''

import numpy as np
from yaff import *
from molmod.io import *
from molmod.units import *
import os

fn_chk = 'opt.chk'
fn_ff = 'pars.txt'
fn_data = 'StrainStress.txt'

conv_gpos = 1e-8
conv_dpos = 1e-6
chk = chk.load_chk(fn_chk)
rvecs = chk['rvecs']
posi = chk['pos']
max_runs = 20000
rcond = 1e-6


def makeMatrix(strain):
    ''' Kan korter, maar hier eenvoudiger te interpreteren '''
    Matrix = np.zeros((6,21))
    Matrix[0,0:6] = strain[:]
    Matrix[1,1] = strain[0]
    Matrix[1,6:11] = strain[1:]
    Matrix[2,2] = strain[0]
    Matrix[2,7] = strain[1]
    Matrix[2,11:15] = strain[2:]
    Matrix[3,3] = strain[0]
    Matrix[3,8] = strain[1]
    Matrix[3,12] = strain[2]
    Matrix[3,15:18] = strain[3:]
    Matrix[4,4] = strain[0]
    Matrix[4,9] = strain[1]
    Matrix[4,13] = strain[2]
    Matrix[4,16] = strain[3]
    Matrix[4,18:20] = strain[4:]
    Matrix[5,5] = strain[0]
    Matrix[5,10] = strain[1]
    Matrix[5,14] = strain[2]
    Matrix[5,17] = strain[3]
    Matrix[5,19] = strain[4]
    Matrix[5,20] = strain[5]
    return Matrix
    
def combineMatrices(strainmatrix):
    m = strainmatrix.shape[1]
    Total = np.zeros((6*m,21))
    for I in xrange(m):
        Total[6*I:6*(I+1),:] = makeMatrix(strainmatrix[:,I])
    return Total

def determineStress(vtens,V):
    stress = np.ones((6))
    stress[0] = 1./V*vtens[0,0]
    stress[1] = 1./V*vtens[1,1]
    stress[2] = 1./V*vtens[2,2]
    stress[3] = 1./V*vtens[1,2]
    stress[4] = 1./V*vtens[0,2]
    stress[5] = 1./V*vtens[0,1]
    return stress
    
def determineCellVector(strain,rvecs):
    rvecs_new = np.ones((3,3))
    strain_matrix = np.zeros((3,3))
    strain_matrix[0,0] = 1. + strain[0]
    strain_matrix[1,1] = 1. + strain[1]
    strain_matrix[2,2] = 1. + strain[2]
    strain_matrix[1,2] = strain[3]/2
    strain_matrix[2,1] = strain[3]/2
    strain_matrix[0,2] = strain[4]/2
    strain_matrix[2,0] = strain[4]/2
    strain_matrix[1,0] = strain[5]/2
    strain_matrix[0,1] = strain[5]/2
    rvecs_new = np.dot(rvecs,strain_matrix)
    return rvecs_new


###### KIEZEN STRAIN ######
''' http://dx.doi.org/10.1063/1.2711762 '''
''' x = +/- 0.007, 0.013 '''

aantalstrains = 48
x1 = 0.003
x2 = 0.006
x3 = 0.009
x4 = 0.012

strainmatrix = np.zeros((6,aantalstrains))

for I in xrange(6):
    strainmatrix[I,I] = x1
    strainmatrix[I,I+6] = -x1
    strainmatrix[I,I+12] = x2
    strainmatrix[I,I+18] = -x2
    strainmatrix[I,I+24] = x3
    strainmatrix[I,I+30] = -x3
    strainmatrix[I,I+36] = x4
    strainmatrix[I,I+42] = -x4
    
A = combineMatrices(strainmatrix)

### Calculating Stress Tensor ###

system = System.from_file(fn_chk)
fns = []
for fn in os.listdir(os.getcwd()):
    if fn.startswith('pars') and fn.endswith('.txt'):
        fns.append(fn)
ff = ForceField.generate(system, fns, rcut=15*angstrom, alpha_scale=3.2, gcut_scale=1.5, smooth_ei=True) 
 
dof = CartesianDOF(ff,gpos_rms=conv_gpos, dpos_rms=conv_dpos) 
opt = CGOptimizer(dof)
opt.run(max_runs)
stress = []
counts = []

for I in xrange(aantalstrains):
    rvecs_temp = determineCellVector(strainmatrix[:,I],rvecs)
    ff.update_rvecs(rvecs_temp)
    ff.update_pos(posi)
    dof = CartesianDOF(ff,gpos_rms=conv_gpos, dpos_rms=conv_dpos) 
    opt = CGOptimizer(dof)
    opt.run(max_runs)
    counts.append(opt.counter)
    ff.compute(vtens=np.ones((3,3)))
    stress.append(determineStress(ff.vtens,system.cell.volume))

y = []

for I in xrange(aantalstrains):
    for J in xrange(6):
        y.append(stress[I][J])
        
Y = np.array(y)

C,res,rank,s = np.linalg.lstsq(A,Y,rcond)

fdata = open(fn_data,'w')

print >> fdata, 'C11 [GPa]:'+'\t'+'\t'+str(C[0]/pascal/10**9)
print >> fdata, 'C12 [GPa]:'+'\t'+'\t'+str(C[1]/pascal/10**9)
print >> fdata, 'C13 [GPa]:'+'\t'+'\t'+str(C[2]/pascal/10**9)
print >> fdata, 'C14 [GPa]:'+'\t'+'\t'+str(C[3]/pascal/10**9)
print >> fdata, 'C15 [GPa]:'+'\t'+'\t'+str(C[4]/pascal/10**9)
print >> fdata, 'C16 [GPa]:'+'\t'+'\t'+str(C[5]/pascal/10**9)
print >> fdata, 'C22 [GPa]:'+'\t'+'\t'+str(C[6]/pascal/10**9)
print >> fdata, 'C23 [GPa]:'+'\t'+'\t'+str(C[7]/pascal/10**9)
print >> fdata, 'C24 [GPa]:'+'\t'+'\t'+str(C[8]/pascal/10**9)
print >> fdata, 'C25 [GPa]:'+'\t'+'\t'+str(C[9]/pascal/10**9)
print >> fdata, 'C26 [GPa]:'+'\t'+'\t'+str(C[10]/pascal/10**9)
print >> fdata, 'C33 [GPa]:'+'\t'+'\t'+str(C[11]/pascal/10**9)
print >> fdata, 'C34 [GPa]:'+'\t'+'\t'+str(C[12]/pascal/10**9)
print >> fdata, 'C35 [GPa]:'+'\t'+'\t'+str(C[13]/pascal/10**9)
print >> fdata, 'C36 [GPa]:'+'\t'+'\t'+str(C[14]/pascal/10**9)
print >> fdata, 'C44 [GPa]:'+'\t'+'\t'+str(C[15]/pascal/10**9)
print >> fdata, 'C45 [GPa]:'+'\t'+'\t'+str(C[16]/pascal/10**9)
print >> fdata, 'C46 [GPa]:'+'\t'+'\t'+str(C[17]/pascal/10**9)
print >> fdata, 'C55 [GPa]:'+'\t'+'\t'+str(C[18]/pascal/10**9)
print >> fdata, 'C56 [GPa]:'+'\t'+'\t'+str(C[19]/pascal/10**9)
print >> fdata, 'C66 [GPa]:'+'\t'+'\t'+str(C[20]/pascal/10**9)
print >> fdata, '--------------------------------------'
print >> fdata, 'Residual (lstsq) = %s'%res
print >> fdata, 'Rank A = %s'%rank
print >> fdata, 'Singuliere waarden A = %s'%s

for I in xrange(aantalstrains):
    print >> fdata, 'Strain %s opt runs: %s'%(I,counts[I])

fdata.close()
