#! /usr/bin/env python

import numpy as np
import h5py
from yaff import System, Cell
from molmod.units import angstrom

f_traj = 'traj_1.h5'
f_out = 'traj_1_strain.h5'
f_clusters = '../clusters.csv'

clusters = np.loadtxt(f_clusters, dtype=int, delimiter=',')

with h5py.File(f_traj, 'r') as f:
    pos = np.array(f['trajectory/pos'])
    rvecs = np.array(f['trajectory/cell'])
    time = np.array(f['trajectory/time'])
    
# array to store the cluster COMs
cluster_pos = np.zeros((len(time),clusters.shape[0],3))
# array to store the virtual clusters
cluster_pos_int = np.zeros((len(time),clusters.shape[0],3))

for t in range(len(time)):
    cell_t = Cell(rvecs[t,:,:])

    # first determine the real clusters
    for i in range(clusters.shape[0]):
        com_t = pos[t,clusters[i,0],:].copy()
        size_cluster = len(clusters[i])
        for j in range(1, size_cluster):
            tmp = pos[t,clusters[i,j],:] - pos[t,clusters[i,0],:]
            cell_t.mic(tmp)
            com_t += tmp/size_cluster
        cluster_pos[t,i,:] = com_t 

    # then add the positions in between those, to generate boxes in which to calculate the strain fields
    # first shift the posiitons of all 'real' clusters by 0.25 times the supercell length in one direction, as a first guess
    cluster_pos_int[t,:,:] = cluster_pos[t,:,:] + 0.25*rvecs[t,0,:]

    for i in range(cluster_pos_int.shape[1]):
        # find the six real clusters closest to the virtual cluster considered here
        distances_all = np.zeros(cluster_pos_int.shape[1])
        for j in range(len(distances_all)):
            tmp = cluster_pos[t,j,:]-cluster_pos_int[t,i,:]
            cell_t.mic(tmp)
            distances_all[j] = np.sqrt((tmp**2).sum())
        # np.argpartition ensures that the k-th element is in sorted position and that all smaller elements are moved before it
        idx_sorted = np.argpartition(distances_all,6)
        com_old = cluster_pos_int[t,i,:].copy()
        for j in idx_sorted[:6]:
            tmp = cluster_pos[t,j,:] - com_old
            cell_t.mic(tmp)
            cluster_pos_int[t,i,:] += tmp/6

# write out the relevant data in HDF5 format

with h5py.File(f_out, 'w') as f:
    fgrp_0 = f.create_group('system')
    fgrp_t = f.create_group('trajectory')
    
    fgrp_0.create_dataset('numbers', data = [36]*cluster_pos.shape[1] + [54]*cluster_pos_int.shape[1])
    fgrp_0.create_dataset('pos', data = np.vstack(( cluster_pos[0,:,:], cluster_pos_int[0,:,:] )))
    fgrp_0.create_dataset('rvecs', data = rvecs[0,:,:])
    
    fgrp_t.create_dataset('pos', data = np.concatenate(( cluster_pos, cluster_pos_int ), axis=1))
    fgrp_t.create_dataset('cell', data = rvecs)
    fgrp_t.create_dataset('time', data=time)
