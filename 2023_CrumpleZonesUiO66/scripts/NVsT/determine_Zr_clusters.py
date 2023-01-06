#! /usr/bin/env python

import numpy as np
from yaff import System, Cell
from molmod.units import angstrom

sys = System.from_file('init.chk')
threshold = 7*angstrom

idx_Zr = np.where(sys.numbers == 40)[0]
pos_Zr = sys.pos[idx_Zr,:]
cell = Cell(sys.cell.rvecs)

n_Zr = len(idx_Zr)
n_cluster = int(n_Zr/6)

mutual_distances = np.zeros((n_Zr,n_Zr))

for i in range(n_Zr):
    pos_i = pos_Zr[i,:]
    mutual_distances[i,i] = 100*angstrom
    for j in range(i+1, n_Zr):
        pos_j = pos_Zr[j,:]
        delta = pos_j-pos_i
        cell.mic(delta)
        mutual_distances[i,j] = np.sqrt(np.sum(delta**2))
        mutual_distances[j,i] = mutual_distances[i,j]

def select_cluster(idx_i, dists, already_chosen):
    cluster_i = [idx_Zr[idx_i]]
    already_chosen.append(idx_i)
    distances = dists[idx_i,:]
    for j, dist in enumerate(distances):
        if dist < threshold: 
            cluster_i.append(idx_Zr[j])
            already_chosen.append(j)
    return cluster_i, already_chosen
    
clusters = []
already_chosen = []
while len(already_chosen) < n_Zr:
    nucleus = 0
    while nucleus in already_chosen:
        nucleus += 1
    cluster_i, already_chosen = select_cluster(nucleus, mutual_distances, already_chosen)
    clusters.append(cluster_i)
    
np.savetxt('clusters.csv', np.array(clusters), delimiter=',', fmt='%d')
