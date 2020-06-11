#! /usr/bin/env python

import h5py
import numpy as np
import matplotlib.pyplot as pt
from molmod.constants import *
from molmod.units import *

def calculate_elastic_tensor(start_time=1*nanosecond, prod_time=10*nanosecond):
    # Retrieve the necessary data
    def rdfile(fname):
        f = h5py.File(fname, 'r')
        time = np.array(f['trajectory/time'])
        volume = np.array(f['trajectory/volume'])
        temp = np.array(f['trajectory/temp'])
        cell = np.array(f['trajectory/cell'])
        press = np.array(f['trajectory/press'])
        f.close()
        return time, volume, temp, cell, press
    time1, volume1, temp1, cell1, press1 = rdfile('traj_1.h5')
    time2, volume2, temp2, cell2, press2 = rdfile('traj_2.h5')
    time3, volume3, temp3, cell3, press3 = rdfile('traj_3.h5')
    time4, volume4, temp4, cell4, press4 = rdfile('traj_4.h5')
    time5, volume5, temp5, cell5, press5 = rdfile('traj_5.h5')
#    time6, volume6, temp6, cell6, press6 = rdfile('traj_6.h5')
#    time7, volume7, temp7, cell7, press7 = rdfile('traj_7.h5')
#    time8, volume8, temp8, cell8, press8 = rdfile('traj_8.h5')
#    time9, volume9, temp9, cell9, press9 = rdfile('traj_9.h5')
#    time10, volume10, temp10, cell10, press10 = rdfile('traj_10.h5')

    time = np.concatenate((time1,time2,time3,time4,time5))##,time6))#,time7,time8,time9,time10))
    volume = np.concatenate((volume1,volume2,volume3,volume4,volume5))#,volume6))#,volume7,volume8,volume9,volume10))
    temp = np.concatenate((temp1,temp2,temp3,temp4,temp5))#,temp6))#,temp7,temp8,temp9,temp10))
    cell = np.concatenate((cell1,cell2,cell3,cell4,cell5))#,cell6))#,cell7,cell8,cell9,cell10))
    press = np.concatenate((press1,press2,press3,press4,press5))#,press6))#,press7,press8,press9,press10))
    eq_steps = int(start_time/(time[1]-time[0]))
    prod_steps = int(prod_time/(time[1]-time[0]))
    if eq_steps + prod_steps > len(time):
	print '%d ns of data requested, but only %.1f available' %((start_time+prod_time)/nanosecond, time[-1]/nanosecond)
        prod_steps = len(time)-eq_steps
    else:
        print 'Sufficient steps found: %.1f ns will be discarded' %((time[-1]-start_time-prod_time)/nanosecond)
    volume = volume[eq_steps:eq_steps+prod_steps]
    temp = temp[eq_steps:eq_steps+prod_steps]
    cell = cell[eq_steps:eq_steps+prod_steps]
    press = press[eq_steps:eq_steps+prod_steps]

    # Calculate the averages and its inverse
    print 'eq = %.1f ns, prod = %.1f ns' %(time[eq_steps]/nanosecond, time[eq_steps+prod_steps-1]/nanosecond)
    cell_av = cell.mean(axis=0)
    temp_av = temp.mean()
    print temp_av
    press_av = press.mean()
    print press_av/(1e9*pascal)
    volume_av = volume.mean()
    cell_av_inv = np.linalg.inv(cell_av)
    print volume_av/angstrom**3

    # Calculate the strain as a function of time
    strain = np.zeros((len(temp),3,3))
    for i in xrange(len(temp)):
        strain[i,:,:] = 0.5*(np.dot(np.dot(cell_av_inv, np.dot(cell[i,:,:], cell[i,:,:].T)), cell_av_inv.T) - np.identity(3))

    # Calculate the full elasticity tensor
    compliance_tensor_full = np.zeros((3,3,3,3))
    for i1 in xrange(3):
        for i2 in xrange(3):
            for i3 in xrange(3):
                for i4 in xrange(3):
                    compliance_tensor_full[i1,i2,i3,i4] = np.cov(strain[:,i1,i2], strain[:,i3,i4])[0,1]
    compliance_tensor_full *= volume_av/(boltzmann*temp_av)

    # Reduce the elasticity tensor to the Voigt notation
    def conv(i):
        if i==0: return 0,0
        if i==1: return 1,1
        if i==2: return 2,2
        if i==3: return 1,2
        if i==4: return 0,2
        if i==5: return 0,1

    elasticity_tensor = np.zeros((6,6))
    compliance_tensor = np.zeros((6,6))
    for i in xrange(6):
        for j in xrange(6):
            i1, i2 = conv(i)
            i3, i4 = conv(j)
            compliance_tensor[i,j] = compliance_tensor_full[i1,i2,i3,i4]

    return np.linalg.inv(compliance_tensor), press_av

def calculate_cubic_constants(tensor):
    c11 = (tensor[0,0] + tensor[1,1] + tensor[2,2])/3.
    c12 = (tensor[0,1] + tensor[0,2] + tensor[1,2])/3.
    c44 = (tensor[3,3] + tensor[4,4] + tensor[5,5])/3.
    return c11, c12, c44

def calculate_orthorhombic_constants(tensor):
    c11 = tensor[0,0]
    c12 = (tensor[0,1] + tensor[1,0])/2.
    c13 = (tensor[0,2] + tensor[2,0])/2.
    c22 = tensor[1,1]
    c23 = (tensor[1,2] + tensor[2,1])/2.
    c33 = tensor[2,2]
    c44 = tensor[3,3]
    c55 = tensor[4,4]
    c66 = tensor[5,5]
    return c11, c12, c13, c22, c23, c33, c44, c55, c66

def print_cubic_born(c11,c12,c44,press):
    print 'C_11 = ' + str(c11/(1e9*pascal)) + ' GPa'
    print 'C_12 = ' + str(c12/(1e9*pascal)) + ' GPa'
    print 'C_44 = ' + str(c44/(1e9*pascal)) + ' GPa'
    print '\n'
    print 'C_11 + 2C_12 + P = ' + str((c11+2*c12+press)/(1e9*pascal)) + ' GPa'
    print 'C_11 - C_12 - 2P = ' + str((c11-c12-2*press)/(1e9*pascal)) + ' GPa'
    print 'C_44 - P = ' + str((c44-press)/(1e9*pascal)) + ' GPa'
    return 0

def print_orthorhombic_born(c11, c12, c13, c22, c23, c33, c44, c55, c66, press):
    print 'C_11 = %s GPa' %str(c11/(1e9*pascal))
    print 'C_12 = %s GPa' %str(c12/(1e9*pascal))
    print 'C_13 = %s GPa' %str(c13/(1e9*pascal))
    print 'C_22 = %s GPa' %str(c22/(1e9*pascal))
    print 'C_23 = %s GPa' %str(c23/(1e9*pascal))
    print 'C_33 = %s GPa' %str(c33/(1e9*pascal))
    print 'C_44 = %s GPa' %str(c44/(1e9*pascal))
    print 'C_55 = %s GPa' %str(c55/(1e9*pascal))
    print 'C_66 = %s GPa' %str(c66/(1e9*pascal))
    print '\n'
    c11p = c11-press
    c22p = c22-press
    c12p = c12+press
    c33p = c33-press
    c13p = c13+press
    c23p = c23+press
    c44p = c44-press
    c55p = c55-press
    c66p = c66-press
    print 'C_11 - P = %s GPa' %str((c11p)/(1e9*pascal))
    print '(C_11 - P)(C_22 - P) - (C_12 + P)**2 = %s GPa**2' %str((c11p*c22p-c12p**2)/(1e9*pascal)**2)
    print '(C_11 - P)(C_22 - P)(C_33 - P) + 2(C_12 + P)(C_13 + P)(C_23 + P) - (C_11 - P)(C_23 + P)**2 - (C_22 - P)(C_13 + P)**2 - (C_33 - P)(C_12 + P)**2 = %s GPa**3' %str((c11p*c22p*c33p+2*c12p*c13p*c23p-c11p*c23p*c23p-c22p*c13p*c13p-c33p*c12p*c12p)/(1e9*pascal)**3)
    print 'C_44 - P = %s GPa' %str(c44p/(1e9*pascal))
    print 'C_55 - P = %s GPa' %str(c55p/(1e9*pascal))
    print 'C_66 - P = %s GPa' %str(c66p/(1e9*pascal))
    return 0

def print_tensor(matrix):
    for i in range(matrix.shape[0]):
        s=''
        for j in range(matrix.shape[1]):
            s += '%.5e\t' %(matrix[i,j])
        print(s)

def print_vector(vector):
    s = ''
    for i in range(len(vector)):
        s += '%.3f\t' %(vector[i])
    print(s)

tensor, press = calculate_elastic_tensor(1*nanosecond, 9*nanosecond)
eigval, eigvec = np.linalg.eig(tensor)
print('\nStiffness tensor C [GPa]')
print_tensor(tensor/(1e9*pascal))
print('\nEigenvalues C [GPa]')
print_vector(eigval/(1e9*pascal))

c11, c12, c44 = calculate_cubic_constants(tensor)
print_cubic_born(c11, c12, c44, press)
