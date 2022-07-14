#!/usr/bin/env python
import os
import h5py as h5
import numpy as np
import matplotlib.pyplot as pt
from molmod.units import angstrom

def determine_CV(fn_traj, indices_mol, indices_plane):
    with h5.File(fn_traj, mode='r') as f:
        numbers_Zn = np.array(f['system/numbers'])[indices_plane]
        assert (numbers_Zn == 30).all()
        pos_plane = np.array(f['trajectory/pos'])[:,indices_plane,:]
        pos_mol = np.array(f['trajectory/pos'])[:,indices_mol,:] # slicing h5 should be in increasing order;

    diffs = np.array([pos_plane[:,(i+2)%6,:]-pos_plane[:,i,:] for i in range(len(indices_plane))]) # array of r13, r24, r35, r46, r51, r62, where rij = r_j - r_i. Note: takes on form 6 x n_time x 3

    # Calculate 6MR normals in the sense of movement of the SF6 molecules
    normals = np.array([np.cross(diffs[i,:,:],diffs[(i+2)%6,:,:]) for i in range(len(indices_plane))]) # construct array of normals r13 x r35, r24 x r46, etc.
    normals_av = np.mean(normals, axis=0)
    normals_av = normals_av/ np.linalg.norm(normals_av,axis=1)[:,np.newaxis]

    # Calculate the collective variable
    r = np.mean(pos_mol, axis=1) - np.mean(pos_plane, axis=1)
    CVs = np.array([np.dot(r[i,:],normals_av[i,:]) for i in range(r.shape[0])])
    
    return CVs


indices_SF6 = np.array([2208,2209,2210,2211,2212,2213,2214])
indices_6MR = np.array([494, 1326, 1324, 778, 768, 500])

bins=np.linspace(-12,12,240)
dirs_25 = []
for i in np.arange(-11,11.5,0.5):
    if i < 0:
        dirs_25.append('../m%d' %(np.around(-10*i)))
    else:
        dirs_25.append('../%d' %(np.around(10*i)))
dirs_50 = []
for i in np.arange(-3,3.5,0.5):
    if i < 0:
        dirs_50.append('../m%d_50' %(np.around(-10*i)))
    else:
        dirs_50.append('../%d_50' %(np.around(10*i)))

dirs_75 = []
for i in np.arange(-3,3.5,0.5):
    if i < 0:
        dirs_75.append('../m%d_75' %(np.around(-10*i)))
    else:
        dirs_75.append('../%d_75' %(np.around(10*i)))

dirs_100 = []
for i in np.arange(-3,3.5,0.5):
    if i < 0:
        dirs_100.append('../m%d_100' %(np.around(-10*i)))
    else:
        dirs_100.append('../%d_100' %(np.around(10*i)))

dirs_150 = []
for i in np.arange(-3,3.5,0.5):
    if i < 0:
        dirs_150.append('../m%d_150' %(np.around(-10*i)))
    else:
        dirs_150.append('../%d_150' %(np.around(10*i)))

dirs_200 = []
for i in np.arange(-3,3.5,0.5):
    if i < 0:
        dirs_200.append('../m%d_200' %(np.around(-10*i)))
    else:
        dirs_200.append('../%d_200' %(np.around(10*i)))

dirs_300 = []
for i in np.arange(-3,3.5,0.5):
    if i < 0:
        dirs_300.append('../m%d_300' %(np.around(-10*i)))
    else:
        dirs_300.append('../%d_300' %(np.around(10*i)))

dirs_400 = []
for i in np.arange(-3,3.5,0.5):
    if i < 0:
        dirs_400.append('../m%d_400' %(np.around(-10*i)))
    else:
        dirs_400.append('../%d_400' %(np.around(10*i)))


fn_name = 'traj_1.h5'
fn_target = 'CVs_1SF6'

if not os.path.isfile('../results/%s.csv'%fn_target):
    CV_all_25 = determine_CV('%s/%s' %(dirs_25[0],fn_name), indices_SF6, indices_6MR)
    for i in range(1,len(dirs_25)):
        print('working on %s' %dirs_25[i])
        CV_all_25 = np.vstack((CV_all_25, determine_CV('%s/%s' %(dirs_25[i],fn_name), indices_SF6, indices_6MR)))
    np.savetxt('../results/%s.csv'%fn_target, CV_all_25, delimiter=',')
else:
    CV_all_25 = np.loadtxt('../results/%s.csv'%fn_target, delimiter=',')

if not os.path.isfile('../results/%s_50.csv'%fn_target):
    CV_all_50 = determine_CV('%s/%s' %(dirs_50[0],fn_name), indices_SF6, indices_6MR)
    for i in range(1,len(dirs_50)):
        print('working on %s' %dirs_50[i])
        CV_all_50 = np.vstack((CV_all_50, determine_CV('%s/%s' %(dirs_50[i],fn_name), indices_SF6, indices_6MR)))
    np.savetxt('../results/%s_50.csv'%fn_target, CV_all_50, delimiter=',')
else:
    CV_all_50 = np.loadtxt('../results/%s_50.csv'%fn_target, delimiter=',')

if not os.path.isfile('../results/%s_75.csv'%fn_target):
    CV_all_75 = determine_CV('%s/%s' %(dirs_75[0],fn_name), indices_SF6, indices_6MR)
    for i in range(1,len(dirs_75)):
        print('working on %s' %dirs_75[i])
        CV_all_75 = np.vstack((CV_all_75, determine_CV('%s/%s' %(dirs_75[i],fn_name), indices_SF6, indices_6MR)))
    np.savetxt('../results/%s_75.csv'%fn_target, CV_all_75, delimiter=',')
else:
    CV_all_75 = np.loadtxt('../results/%s_75.csv'%fn_target, delimiter=',')

if not os.path.isfile('../results/%s_100.csv'%fn_target):
    CV_all_100 = determine_CV('%s/%s' %(dirs_100[0],fn_name), indices_SF6, indices_6MR)
    for i in range(1,len(dirs_100)):
        print('working on %s' %dirs_100[i])
        CV_all_100 = np.vstack((CV_all_100, determine_CV('%s/%s' %(dirs_100[i],fn_name), indices_SF6, indices_6MR)))
    np.savetxt('../results/%s_100.csv'%fn_target, CV_all_100, delimiter=',')
else:
    CV_all_100 = np.loadtxt('../results/%s_100.csv'%fn_target, delimiter=',')

if not os.path.isfile('../results/%s_150.csv'%fn_target):
    CV_all_150 = determine_CV('%s/%s' %(dirs_150[0],fn_name), indices_SF6, indices_6MR)
    for i in range(1,len(dirs_150)):
        print('working on %s' %dirs_150[i])
        CV_all_150 = np.vstack((CV_all_150, determine_CV('%s/%s' %(dirs_150[i],fn_name), indices_SF6, indices_6MR)))
    np.savetxt('../results/%s_150.csv'%fn_target, CV_all_150, delimiter=',')
else:
    CV_all_150 = np.loadtxt('../results/%s_150.csv'%fn_target, delimiter=',')

if not os.path.isfile('../results/%s_200.csv'%fn_target):
    CV_all_200 = determine_CV('%s/%s' %(dirs_200[0],fn_name), indices_SF6, indices_6MR)
    for i in range(1,len(dirs_200)):
        print('working on %s' %dirs_200[i])
        CV_all_200 = np.vstack((CV_all_200, determine_CV('%s/%s' %(dirs_200[i],fn_name), indices_SF6, indices_6MR)))
    np.savetxt('../results/%s_200.csv'%fn_target, CV_all_200, delimiter=',')
else:
    CV_all_200 = np.loadtxt('../results/%s_200.csv'%fn_target, delimiter=',')

if not os.path.isfile('../results/%s_300.csv'%fn_target):
    CV_all_300 = determine_CV('%s/%s' %(dirs_300[0],fn_name), indices_SF6, indices_6MR)
    for i in range(1,len(dirs_300)):
        print('working on %s' %dirs_300[i])
        CV_all_300 = np.vstack((CV_all_300, determine_CV('%s/%s' %(dirs_300[i],fn_name), indices_SF6, indices_6MR)))
    np.savetxt('../results/%s_300.csv'%fn_target, CV_all_300, delimiter=',')
else:
    CV_all_300 = np.loadtxt('../results/%s_300.csv'%fn_target, delimiter=',')

if not os.path.isfile('../results/%s_400.csv'%fn_target):
    CV_all_400 = determine_CV('%s/%s' %(dirs_400[0],fn_name), indices_SF6, indices_6MR)
    for i in range(1,len(dirs_400)):
        print('working on %s' %dirs_400[i])
        CV_all_400 = np.vstack((CV_all_400, determine_CV('%s/%s' %(dirs_400[i],fn_name), indices_SF6, indices_6MR)))
    np.savetxt('../results/%s_400.csv'%fn_target, CV_all_400, delimiter=',')
else:
    CV_all_400 = np.loadtxt('../results/%s_400.csv'%fn_target, delimiter=',')

CV_all = np.vstack((CV_all_25, CV_all_50, CV_all_75, CV_all_100, CV_all_150, CV_all_200, CV_all_300, CV_all_400))
 
pt.clf()
pt.xlabel('CV (angstrom)')
pt.ylabel('probability')
for i in range(CV_all.shape[0]):
    hist, bin_edges = np.histogram(CV_all[i,:]/angstrom, bins=bins, density=True)
    pt.plot( 0.5*(bin_edges[1:] + bin_edges[:-1]), hist)
pt.savefig('../results/%s.svg' %fn_target)
