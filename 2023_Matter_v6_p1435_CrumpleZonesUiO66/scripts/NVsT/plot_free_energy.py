#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as pt
from molmod.units import *

def press_fit(volume_exp, polycoef, scale):
    fit_press = np.zeros(len(volume_exp))
    for i in range(len(polycoef)):
        fit_press += polycoef[len(polycoef)-i-1]*(volume_exp/scale)**i
    return fit_press     

def bulk_modulus(volume_opt, polycoef, scale):
    bulk_mod = 0
    for i in range(len(polycoef)):
        bulk_mod -=  i*polycoef[len(polycoef)-i-1]*(volume_opt/scale)**i
    return bulk_mod

def fit_data(fname, n_unitcell=1):
    data = np.genfromtxt(fname, delimiter = '\t')
    volume = data[:,0]*angstrom**3/n_unitcell
    press = data[:,1]*(1e6*pascal)
    h = np.argmin(np.abs(volume-6000*angstrom**3))
    volume = volume[h:]
    press = press[h:]
    poly = np.polyfit(volume/np.mean(volume), press, 7)
    fit = press_fit(volume, poly, np.mean(volume))

    return volume, press, poly, fit

def trans_press(fit):
    h = len(fit)/2
    max_press = np.amax(fit[h:])
    min_press = np.amin(fit[:h])
    return np.array([min_press, max_press])/(1e6*pascal)

def extr_vols(fit, vols, poly):
    def index_of_abs_min(fit_sub):
        abs_min = min(fit_sub, key=abs)
        index = np.where(fit_sub == abs_min)[0]
        return index[0]

    vol_opt = vols[index_of_abs_min(fit)]

    bulk_opt = bulk_modulus(vol_opt, poly, np.mean(vols))
    return vol_opt/angstrom**3, bulk_opt/(1e9*pascal)

def free_ener_fit(volume_exp, polycoef, scale):
    fit_free_ener = np.zeros(len(volume_exp))
    for i in range(len(polycoef)):
        fit_free_ener -= polycoef[len(polycoef)-i-1]*(volume_exp/scale)**i*volume_exp/(i+1)
    return fit_free_ener - np.amin(fit_free_ener)  

def determine_loc(fit, volume):
    h = int(len(fit)/2)
    return volume[h+np.argmax(fit[h:])]/angstrom**3, np.amax(fit[h:])/(1e6*pascal)
	
	
n_unitcell=2*2*2
vol, press, poly, fit = fit_data('./PvsV_eq100ps_av900ps.csv', n_unitcell)
volume_finegrid = np.arange(np.amin(vol), np.amax(vol), 0.1*angstrom**3)
fit = press_fit(volume_finegrid, poly, np.mean(vol))
ener = free_ener_fit(volume_finegrid, poly, np.mean(vol))
opt_vol, bulk = extr_vols(fit, volume_finegrid, poly)
loc_volume, loc_pressure = determine_loc(press, vol)

pt.clf()
pt.figure(figsize=(8, 3))
pt.xlabel(r'Cell volume [$\rm{\AA}^3$]')
pt.ylabel(r'Pressure [MPa]')
pt.title(r'Pressure profile of UiO-CN12 at 300 K')
pt.xlim([6000,9500])
pt.ylim([-1000,2000])

pt.grid(b=True, which = 'major', color='0.25', linestyle=':')

pt.plot((vol/angstrom**3), press/(1e6*pascal), linewidth=0, marker='.', markersize=10, color='#3690c0')
xmin, xmax = pt.xlim()
ymin, ymax = pt.ylim()
y_range = ymax-ymin
x_range = xmax-xmin

# Loss of crystallinity information
pt.plot((loc_volume - 0.1*x_range, loc_volume + 0.1*x_range), (loc_pressure, loc_pressure), linewidth=1, color='#e31a1c')
pt.plot((loc_volume, loc_volume), (loc_pressure - 0.05*y_range, loc_pressure + 0.05*y_range), linewidth=1, color='#e31a1c')
pt.text(loc_volume + 0.02*x_range, loc_pressure+0.06*y_range, r'%.0f $\rm{\AA}^3$ (%.2f $\rm{\AA}$)'  %(loc_volume, (loc_volume)**(1./3)), color='#e31a1c', bbox={'facecolor':'none', 'linewidth':0})
pt.text(loc_volume + 0.02*x_range, loc_pressure+0.02*y_range, r'%.0f MPa'  %loc_pressure, color='#e31a1c', bbox={'facecolor':'none', 'linewidth':0})

# Equilibrium information
pt.plot((opt_vol - 0.1*x_range, opt_vol + 0.1*x_range), (0, 0), linewidth=1, color='#238b45')
pt.plot((opt_vol, opt_vol), (-0.05*y_range, 0.05*y_range), linewidth=1, color='#238b45')
pt.text(opt_vol + 0.02*x_range, 0.06*y_range, r'%.0f $\rm{\AA}^3$ (%.2f $\rm{\AA}$)'  %(opt_vol, (opt_vol)**(1./3)), color='#238b45', bbox={'facecolor':'none', 'linewidth':0})
pt.text(opt_vol + 0.02*x_range, 0.02*y_range, r'$K = $ %.1f GPa' %bulk, color='#238b45', bbox={'facecolor':'none', 'linewidth':0})

#pt.plot((volume_finegrid/angstrom**3), fit/(1e6*pascal), color='#034e7b')

pt.savefig('PvsV_UiO_CN12_all.pdf', format='pdf', bbox_inches = 'tight')
pt.savefig('PvsV_UiO_CN12_all.svg', bbox_inches = 'tight')
