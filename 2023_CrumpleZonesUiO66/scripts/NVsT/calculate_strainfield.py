#! /usr/bin/env python

import numpy as np
import h5py
from yaff import System, Cell
from molmod.units import angstrom, picosecond
import matplotlib.pyplot as pt
import matplotlib.animation

def get_correct_order(f_index, key='trajectory/sorting_indices'):
    # This function orders all nodes in such a way that the nodes 
    # are ordered  in increasing order of x, y, and z
    
    with h5py.File(f_index, 'a') as h:
        if key in h:
            # simply take the ordering calculated before
            ind_all = np.array(h[key]).astype(int)
            print('Reusing previously calculated ordering')
        else:
            print('Calculating correct ordering')
            # calculate the correct ordering
            pos_0 = np.array(h['trajectory/pos'][0,:,:]).copy()
            n_subcells_1dim = int(round(pos_0.shape[0]**(1./3)))
            rvecs_0 = np.array(h['trajectory/cell'][0,:,:])
            
            # first get the original positions within mic
            cell_0 = Cell(rvecs_0)
            for i in range(pos_0.shape[0]):
                tmp = pos_0[i,:]
                cell_0.mic(tmp)
                for j in range(3):
                    if tmp[j] < 2*angstrom: tmp += cell_0.rvecs[j,:]
                pos_0[i,:] = tmp
                
            # order sequentially over z, y, and x (sequential order ensures similar but unequal values are binned together)
            ind_1 = np.argsort(pos_0[:,2])
            pos_0_sorted = pos_0[ind_1,:]
            tmp = [np.argsort(pos_0_sorted[n_subcells_1dim**2*i:n_subcells_1dim**2*(i+1),1])+i*n_subcells_1dim**2 for i in range(n_subcells_1dim)]
            ind_2 = np.concatenate([tmp[i] for i in range(n_subcells_1dim)])
            pos_0_sorted = pos_0_sorted[ind_2,:]
            tmp = [np.argsort(pos_0_sorted[n_subcells_1dim*i:n_subcells_1dim*(i+1),0])+i*n_subcells_1dim for i in range(n_subcells_1dim**2)]
            ind_3 = np.concatenate([tmp[i] for i in range(n_subcells_1dim**2)])
            pos_0_sorted = pos_0_sorted[ind_3,:]
            ind_all = ind_1[ind_2[ind_3]]
            
            # also store sorting indices for subsequent use
            hgrp_t = h['trajectory']
            hgrp_t.create_dataset('sorting_indices', data = ind_all)

    return ind_all

def index_to_100(i, n_subcells_1dim):
    return (i+1)%(n_subcells_1dim)-i%(n_subcells_1dim)

def index_to_010(i, n_subcells_1dim):
    return (i+n_subcells_1dim)%(n_subcells_1dim**2)-i%(n_subcells_1dim**2)

def index_to_001(i, n_subcells_1dim):
    return (i+n_subcells_1dim**2)%(n_subcells_1dim**3)-i%(n_subcells_1dim**3)
    
def create_subcells(pos_sorted, rvecs):
    # This function creates subcells from the given node positions,
    # returning their (average) rvecs and center-of-mass position
    # as n_time x n_subcells x 3 x 3 and n_time x n_subcells x 3
    # NumPy arrays, respectively
    
    # initialize the output arrays
    subcells_rvecs = np.zeros((pos_sorted.shape[0],pos_sorted.shape[1],3,3))
    subcells_com = np.zeros((pos_sorted.shape[0],pos_sorted.shape[1],3))
    
    # determine number of subcells in one direction (assuming same number of nodes in all three directions)
    n_subcells_1dim = int(round(pos_sorted.shape[1]**(1./3)))
    
    # first create the relative cell vectors, without mic at this point, for all time steps at the same time
    for i in range(pos_sorted.shape[1]): # run over each node, for all timesteps at the same time
        pos_000 = pos_sorted[:,i,:]
        pos_100 = pos_sorted[:,i+index_to_100(i,n_subcells_1dim),:]
        pos_010 = pos_sorted[:,i+index_to_010(i,n_subcells_1dim),:]
        pos_110 = pos_sorted[:,i+index_to_100(i,n_subcells_1dim)+index_to_010(i,n_subcells_1dim),:]
        pos_001 = pos_sorted[:,i+index_to_001(i,n_subcells_1dim),:]
        pos_101 = pos_sorted[:,i+index_to_100(i,n_subcells_1dim)+index_to_001(i,n_subcells_1dim),:]
        pos_011 = pos_sorted[:,i+index_to_010(i,n_subcells_1dim)+index_to_001(i,n_subcells_1dim),:]
        pos_111 = pos_sorted[:,i+index_to_100(i,n_subcells_1dim)+index_to_010(i,n_subcells_1dim)+index_to_001(i,n_subcells_1dim),:]

        xvecs_all = np.array([pos_100-pos_000, pos_101-pos_001, pos_110-pos_010, pos_111-pos_011])
        yvecs_all = np.array([pos_010-pos_000, pos_110-pos_100, pos_011-pos_001, pos_111-pos_101])
        zvecs_all = np.array([pos_001-pos_000, pos_101-pos_100, pos_011-pos_010, pos_111-pos_110])

        for t in range(pos_sorted.shape[0]):
            cell_t = Cell(rvecs[t,:,:])
            for j in range(4): # take average over four equivalent vectors in any given direction
                tmp = xvecs_all[j,t,:]
                cell_t.mic(tmp)
                subcells_rvecs[t,i,0,:] += tmp/4

                tmp = yvecs_all[j,t,:]
                cell_t.mic(tmp)
                subcells_rvecs[t,i,1,:] += tmp/4
                
                tmp = zvecs_all[j,t,:]
                cell_t.mic(tmp)
                subcells_rvecs[t,i,2,:] += tmp/4         

            subcells_com[t,i,:] = pos_000[t,:]
            for pos_tmp in [pos_100[t,:], pos_010[t,:], pos_001[t,:], pos_110[t,:], pos_101[t,:], pos_011[t,:], pos_111[t,:]]:
                tmp = pos_tmp - pos_000[t,:]
                cell_t.mic(tmp)
                subcells_com[t,i,:] += tmp/8

    # make a cubic-like visualization by using mic
    mic_move = np.zeros((pos_sorted.shape[1],3))
    
    cell_0 = Cell(rvecs[0,:,:])
    rvecs_invT = np.linalg.inv(rvecs[0,:,:]).T
    already_covered = [0] # indices of the nodes that have been treated; starts from the node at index 0
    
    subcells_com_0 = subcells_com[0,:,:].copy()
    
    for iz in range(n_subcells_1dim):
        for iy in range(n_subcells_1dim):
            for ix in range(n_subcells_1dim):
                i = iz*n_subcells_1dim**2 + iy*n_subcells_1dim + ix
                # each node defines three other nodes, namely in 100, 010, and 001 directions
                for i_next in [i+index_to_100(i, n_subcells_1dim), i+index_to_010(i, n_subcells_1dim), i+index_to_001(i, n_subcells_1dim)]:
                    if i_next in already_covered: pass
                    else:
                        pos_old = subcells_com_0[i_next,:] - subcells_com_0[i,:]
                        pos_new = pos_old.copy()
                        cell_0.mic(pos_new)
                        mic_move[i_next,:] = np.around(np.dot(rvecs_invT, pos_new-pos_old))
                        already_covered.append(i_next)
                        subcells_com_0[i_next,:] += pos_new - pos_old

    mic_move = mic_move.astype(int)
    
    subcells_com2 = subcells_com.copy()
    
    for t in range(subcells_com.shape[0]):
        subcells_com[t,:,:] += np.dot(mic_move, rvecs[t,:,:])
    
    for i in range(subcells_com.shape[1]):
        if np.sum(np.abs(mic_move[i,:])) > 0.5:
            print('moved (%.2f, %.2f, %.2f)\tto\t(%.2f, %.2f, %.2f) \t along (%d, %d, %d)' %(subcells_com2[0,i,0]/angstrom, subcells_com2[0,i,1]/angstrom,subcells_com2[0,i,2]/angstrom, subcells_com[0,i,0]/angstrom, subcells_com[0,i,1]/angstrom,subcells_com[0,i,2]/angstrom, mic_move[i,0], mic_move[i,1], mic_move[i,2]))
        else:
            print('kept (%.2f, %.2f, %.2f)' %(subcells_com2[0,i,0]/angstrom, subcells_com2[0,i,1]/angstrom,subcells_com2[0,i,2]/angstrom))
         
    return subcells_rvecs, subcells_com

def calc_reference_cell(f):
    # calculates the average subcell matrix, obtained by averaging over
    # each timestep and subcell
    subcells_rvecs_ref = np.array(f['trajectory/subcells_rvecs'])
    return np.average(subcells_rvecs_ref, axis=(0,1))
    
def calc_strain(h, h_ref):
    # function to calculate the Lagrangian strain starting from a time series of subcell matrices h
    # (size n_time x n_subcells x 3 x 3) and a reference cell matrix (size 3 x 3)
    
    # calculate h h^T for each timestep and subcell at the same time
    hT = np.swapaxes(h,2,3)
    hhT = np.matmul(h,hT)
    
    # calculate the strain
    h_refinv = np.linalg.inv(h_ref)
    strain = 0.5*(np.matmul(h_refinv, np.matmul(hhT, h_refinv.T)) - np.tile(np.eye(3), (h.shape[0],h.shape[1],1,1)))
    return strain
    
def get_subcell_visualization(subcell_pos, subcell_com_0, rvecs):
    # returns all subcell positions surronding each of the subcells COMs,
    # which implies that the positions at the edge will be visualized twice
    # also returns the pairs of indices along that should be connected through lines
    
    cell_0 = Cell(rvecs[0,:,:])
    rvecs_invT = np.linalg.inv(rvecs[0,:,:]).T
    n_subcells_1dim = int(round(subcell_pos.shape[1]**(1./3)))
    subcell_nodes = np.zeros((subcell_pos.shape[0], (n_subcells_1dim+1)**3, 3))
    subcell_nodes[:] = np.nan
    subcell_lines = np.zeros((3*(n_subcells_1dim)**3+3*2*(n_subcells_1dim)**2+3*n_subcells_1dim, 2), dtype=int) # will contain the node *indices*
    
    counter_nodes = 0
    counter_lines = 0
    
    threshold = 0.1*angstrom
    
    # run over each subcell COM to define the positions of the eight surrounding nodes
    # only the situation at t=0 is necessary for reference
    
    for i in range(subcell_com_0.shape[0]):
        # for the given subcell, get the nodes surrouding it. With the sorted positions,
        # the subcell index corresponds to the index of its node with minimum x, y, z value
        pos_000 = pos_sorted[:,i,:]
        pos_100 = pos_sorted[:,i+index_to_100(i,n_subcells_1dim),:]
        pos_010 = pos_sorted[:,i+index_to_010(i,n_subcells_1dim),:]
        pos_110 = pos_sorted[:,i+index_to_100(i,n_subcells_1dim)+index_to_010(i,n_subcells_1dim),:]
        pos_001 = pos_sorted[:,i+index_to_001(i,n_subcells_1dim),:]
        pos_101 = pos_sorted[:,i+index_to_100(i,n_subcells_1dim)+index_to_001(i,n_subcells_1dim),:]
        pos_011 = pos_sorted[:,i+index_to_010(i,n_subcells_1dim)+index_to_001(i,n_subcells_1dim),:]
        pos_111 = pos_sorted[:,i+index_to_100(i,n_subcells_1dim)+index_to_010(i,n_subcells_1dim)+index_to_001(i,n_subcells_1dim),:]
        
        surrounding_pos = np.array([pos_000, pos_100, pos_010, pos_110, pos_001, pos_101, pos_011, pos_111])
        surrounding_indices = np.array([np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan])
        
        # for each of these nodes, we need to determine whether they actually surround the COM,
        # or whether they were moved because of periodic image conditions
        # in case of the latter, they need to be moved back, and this for every time instant
        
        for idx_pos, pos_old in enumerate(surrounding_pos):
            diff_old = pos_old[0,:] - subcell_com_0[i,:] # calculate original distance between node and COM
            diff_new = diff_old.copy()
            cell_0.mic(diff_new)       
            mic_move = np.around(np.dot(rvecs_invT, diff_new-diff_old))
            
            pos_new = pos_old + np.einsum('i,tik->tk', mic_move, rvecs)
            
            # check whether the position (with possible translation along the cell vectors) is already encountered
            # and keep track of the index where the position is stored to draw the lines afterwards
            
            dist_min = 100*angstrom
            idx_to_store = np.nan
            
            for idx_pos_0, pos_0 in enumerate(subcell_nodes[0,:,:]):
                if ((pos_new[0]-pos_0)**2).sum() < dist_min**2:
                    dist_min = np.sqrt(((pos_new[0]-pos_0)**2).sum())
                    idx_to_store = idx_pos_0
            if dist_min >= threshold:
                subcell_nodes[:,counter_nodes,:] = pos_new.copy()
                idx_to_store = counter_nodes
                counter_nodes += 1
            surrounding_indices[idx_pos] = int(idx_to_store)
            
        # make a list of all twelve potential lines that can be drawn between two adjacent nodes in the box
        # based on the ordering in surrounding_pos (and hence surrounding_indices):
        # 000-100, 000-010, 000-001, 100-110, 100-101, 010-110, 010-011, 001-101, 001-011, 011-111, 101-111, 110-111
        all_lines = np.array([[surrounding_indices[0],surrounding_indices[1]], [surrounding_indices[0],surrounding_indices[2]], [surrounding_indices[0],surrounding_indices[4]], [surrounding_indices[1],surrounding_indices[3]], [surrounding_indices[1],surrounding_indices[5]], [surrounding_indices[2],surrounding_indices[3]], [surrounding_indices[2],surrounding_indices[6]], [surrounding_indices[4],surrounding_indices[5]], [surrounding_indices[4],surrounding_indices[6]], [surrounding_indices[6],surrounding_indices[7]], [surrounding_indices[5],surrounding_indices[7]], [surrounding_indices[3],surrounding_indices[7]]])
        
        for idx_line, line in enumerate(all_lines):
            already_present = False
            for idx_line_0, line_0 in enumerate(subcell_lines):
                if (line[0] == line_0[0] and line[1] == line_0[1]) or (line[0] == line_0[1] and line[1] == line_0[0]):
                    # connection already present
                    already_present = True
                    break
            if not already_present:
                subcell_lines[counter_lines,:] = line
                counter_lines += 1
            
    return subcell_nodes, subcell_lines

# input that is, in principle, independent of the current volume
f_index = 'traj_9275_strain.h5' # file from which the sorting of the nodes is determined
f_refcell = 'traj_9275_strain.h5' # file from which the reference cell is defined

# volume-dependent input and output
f_strain = 'traj_9275_strain.h5' # file for which the strain field is calculated and visualized
f_output = 'traj_9275_strain.mp4' # output animation
animation = True # whether or not to plot a video
f_figname = 'traj_9275' # prefix of the output figures

# load the relevant data from the HDF5 file
with h5py.File(f_strain, 'r') as f:
    pos = np.array(f['trajectory/pos'])
    rvecs = np.array(f['trajectory/cell'])
    time = np.array(f['trajectory/time'])

# order the node positions
ind_all = get_correct_order(f_index)
pos_sorted = pos[:,ind_all,:]

# create subcells, for which the cell vectors are stored as rows, unless they have been created already
with h5py.File(f_strain, 'a') as f:
    if 'trajectory/subcells_rvecs' in f and 'trajectory/subcells_com' in f:
        subcells_rvecs = np.array(f['trajectory/subcells_rvecs'])
        subcells_com = np.array(f['trajectory/subcells_com'])
    else:
        subcells_rvecs, subcells_com = create_subcells(pos_sorted, rvecs)       
        fgrp_t = f['trajectory']
        fgrp_t.create_dataset('subcells_rvecs', data = subcells_rvecs)
        fgrp_t.create_dataset('subcells_com', data = subcells_com)
        
# calculate reference cell, unless it has been calculated already
with h5py.File(f_refcell, 'a') as f:
    if 'trajectory/cell_ref' in f:
        h_ref = np.array(f['trajectory/cell_ref'])
    else:
        h_ref = calc_reference_cell(f)
        fgrp_t = f['trajectory']
        fgrp_t.create_dataset('cell_ref', data = h_ref)
    
# calculate strain field, unless it has been calculated already
with h5py.File(f_strain, 'a') as f:
    if 'trajectory/subcells_strain' in f:
        strain = np.array(f['trajectory/subcells_strain'])
    else:
        strain = calc_strain(subcells_rvecs, h_ref)
        fgrp_t = f['trajectory']
        fgrp_t.create_dataset('subcells_strain', data = strain)
        
# calculate positions of all cell nodes used to calculate strain (with duplicates if necessary),
# including connecting lines
subcell_nodes, subcell_lines = get_subcell_visualization(pos_sorted, subcells_com[0,:,:], rvecs)    

def update_graph(num):
    graph1._offsets3d=(subcell_nodes[num,:,0]/angstrom, subcell_nodes[num,:,1]/angstrom, subcell_nodes[num,:,2]/angstrom)
    graph2pos._offsets3d=(subcells_com[num,mask[num,:],0]/angstrom, subcells_com[num,mask[num,:],1]/angstrom, subcells_com[num,mask[num,:],2]/angstrom)
    graph2neg._offsets3d=(subcells_com[num,np.invert(mask[num,:]),0]/angstrom, subcells_com[num,np.invert(mask[num,:]),1]/angstrom, subcells_com[num,np.invert(mask[num,:]),2]/angstrom)
    graph2pos.set_array(c_value_pos[num,mask[num,:]])
    graph2neg.set_array(c_value_pos[num,np.invert(mask[num,:])])
    for i, line in enumerate(graph_lines):
        line[0].set_data_3d( [subcell_nodes[num,subcell_lines[i,0],0]/angstrom,subcell_nodes[num,subcell_lines[i,1],0]/angstrom], [subcell_nodes[num,subcell_lines[i,0],1]/angstrom,subcell_nodes[num,subcell_lines[i,1],1]/angstrom], [subcell_nodes[num,subcell_lines[i,0],2]/angstrom,subcell_nodes[num,subcell_lines[i,1],2]/angstrom] )
    title.set_text('Strain field, time=%.2f ps' %(time[num]/picosecond))
    
def get_scalar_strain(strain):
    # input: t x n_subcell x 3 x 3 strain matrix containing
    # the strain of each subcell at each timestep
    # output: t x n_subcell matrix containing a linear variant of the strain
    
    #lin_strain = strain[:,:,1,2] # returns yz component of the strain
    lin_strain = strain[:,:,1,1] # returns yy component of the strain
    return lin_strain




if animation:
    fig = pt.figure()
    
    # figure settings
    x_max = 50 # in angstrom
    y_max = 50 # in angstrom
    z_max = 50 # in angstrom
    cmin = -2
    cmax = 0

    ax = pt.axes(projection='3d')
    ax.set_xlim([0,x_max])
    ax.set_xlabel('x (angstrom)')
    ax.set_ylim([0,y_max])
    ax.set_ylabel('y (angstrom)')
    ax.set_zlim([0,z_max])
    ax.set_zlabel('z (angstrom)')
    title = ax.set_title('Strain field')
    ax.view_init(0,0)
    ax.set_box_aspect((1,1,1))

    lin_strain = get_scalar_strain(strain)
    c_value_pos = np.log10(np.abs(lin_strain))
    mask = np.sign(lin_strain) > 0

    graph1 = ax.scatter(subcell_nodes[0,:,0]/angstrom, subcell_nodes[0,:,1]/angstrom, subcell_nodes[0,:,2]/angstrom, color='0.3')
    graph2pos = ax.scatter(subcells_com[0,mask[0,:],0]/angstrom, subcells_com[0,mask[0,:],1]/angstrom, subcells_com[0,mask[0,:],2]/angstrom, s=60, c = c_value_pos[0,mask[0,:]], cmap='YlOrRd')
    graph2neg = ax.scatter(subcells_com[0,np.invert(mask[0,:]),0]/angstrom, subcells_com[0,np.invert(mask[0,:]),1]/angstrom, subcells_com[0,np.invert(mask[0,:]),2]/angstrom, s=60, c = c_value_pos[0,np.invert(mask[0,:])], cmap='YlGnBu')
    graph_lines = [ ax.plot( [subcell_nodes[0,subcell_lines[i,0],0]/angstrom,subcell_nodes[0,subcell_lines[i,1],0]/angstrom], [subcell_nodes[0,subcell_lines[i,0],1]/angstrom,subcell_nodes[0,subcell_lines[i,1],1]/angstrom], [subcell_nodes[0,subcell_lines[i,0],2]/angstrom,subcell_nodes[0,subcell_lines[i,1],2]/angstrom], '--', color='0.3' ) for i in range(subcell_lines.shape[0])]

    graph2pos.set_clim(cmin,cmax)
    graph2neg.set_clim(cmin,cmax)
    cbpos = pt.colorbar(graph2pos)
    cbneg = pt.colorbar(graph2neg)
    cbpos.ax.set_title('pos')
    cbneg.ax.set_title('neg')

    ticks_pos = np.linspace(cmin, cmax, 5, endpoint=True)
    ticks_neg = np.linspace(cmin, cmax, 5, endpoint=True)

    cbpos.set_ticks(ticks_pos)
    cbneg.set_ticks(ticks_neg)
    cbneg.ax.invert_yaxis()

    ani = matplotlib.animation.FuncAnimation(fig, update_graph, pos.shape[0])
    #pt.show()
    ani.save(f_output, fps = 15, extra_args = ['-vcodec', 'libx264'])
    pt.savefig(f_figname + '.png')


fig = pt.figure()

# figure of average strain
x_max = 50 # in angstrom
y_max = 50 # in angstrom
z_max = 50 # in angstrom
cmin = -2
cmax = 0

ax = pt.axes(projection='3d')
ax.set_xlim([0,x_max])
ax.set_xlabel('x (angstrom)')
ax.set_ylim([0,y_max])
ax.set_ylabel('y (angstrom)')
ax.set_zlim([0,z_max])
ax.set_zlabel('z (angstrom)')
title = ax.set_title('Strain field')
ax.view_init(0,0)
ax.set_box_aspect((1,1,1))


lin_strain = get_scalar_strain(strain)
c_value_pos = np.log10(np.abs(lin_strain))
mask = np.sign(lin_strain) > 0

graph1 = ax.scatter((subcell_nodes[:,:,0].mean(axis=0))/angstrom, (subcell_nodes[:,:,1].mean(axis=0))/angstrom, (subcell_nodes[:,:,2].mean(axis=0))/angstrom, color='0.3')
graph2pos = ax.scatter((subcells_com[:,mask[0,:],0].mean(axis=0))/angstrom, (subcells_com[:,mask[0,:],1].mean(axis=0))/angstrom, (subcells_com[:,mask[0,:],2].mean(axis=0))/angstrom, s=60, c = c_value_pos[:,mask[0,:]].mean(axis=0), cmap='YlOrRd')
graph2neg = ax.scatter((subcells_com[:,np.invert(mask[0,:]),0].mean(axis=0))/angstrom, (subcells_com[:,np.invert(mask[0,:]),1].mean(axis=0))/angstrom, (subcells_com[:,np.invert(mask[0,:]),2].mean(axis=0))/angstrom, s=60, c = c_value_pos[:,np.invert(mask[0,:])].mean(axis=0), cmap='YlGnBu')
graph_lines = [ ax.plot( [(subcell_nodes[:,subcell_lines[i,0],0].mean(axis=0))/angstrom,(subcell_nodes[:,subcell_lines[i,1],0].mean(axis=0))/angstrom], [(subcell_nodes[:,subcell_lines[i,0],1].mean(axis=0))/angstrom,(subcell_nodes[:,subcell_lines[i,1],1].mean(axis=0))/angstrom], [(subcell_nodes[:,subcell_lines[i,0],2].mean(axis=0))/angstrom,(subcell_nodes[:,subcell_lines[i,1],2].mean(axis=0))/angstrom], '--', color='0.3' ) for i in range(subcell_lines.shape[0])]

graph2pos.set_clim(cmin,cmax)
graph2neg.set_clim(cmin,cmax)
cbpos = pt.colorbar(graph2pos)
cbneg = pt.colorbar(graph2neg)
cbpos.ax.set_title('pos')
cbneg.ax.set_title('neg')

ticks_pos = np.linspace(cmin, cmax, 5, endpoint=True)
ticks_neg = np.linspace(cmin, cmax, 5, endpoint=True)

cbpos.set_ticks(ticks_pos)
cbneg.set_ticks(ticks_neg)
cbneg.ax.invert_yaxis()

pt.savefig(f_figname + '_av_x.svg')

ax.view_init(90,-90)
pt.savefig(f_figname + '_av_z.svg')
