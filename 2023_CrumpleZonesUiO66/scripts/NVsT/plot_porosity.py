#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as pt

all_data = np.loadtxt('porosity.csv', delimiter=',')

uc_vol = all_data[:,0]/8
density = all_data[:,1]
av_a3 = all_data[:,2]
av_volfrac = all_data[:,3]
av_cm3_g = all_data[:,4]

pt.clf()
pt.figure(figsize=(8, 2))
pt.xlabel(r'Cell volume [$\rm{\AA}^3$]')
pt.ylabel(r'Accessible volume [$\rm{\AA}^3$]')
pt.xlim([6000,9500])
pt.plot(uc_vol, av_a3, linewidth=0, marker='.', markersize=10, color='#31a354')
pt.savefig('porosity_av_a3.pdf', format='pdf', bbox_inches = 'tight')
pt.savefig('porosity_av_a3.svg', bbox_inches = 'tight')

pt.clf()
pt.figure(figsize=(8, 2))
pt.xlabel(r'Cell volume [$\rm{\AA}^3$]')
pt.ylabel(r'Accessible volume [\%]')
pt.xlim([6000,9500])
pt.ylim([0,0.3])
pt.plot(uc_vol[av_volfrac>0.1], av_volfrac[av_volfrac>0.1], linewidth=0, marker='.', markersize=10, color='#31a354')
pt.plot(uc_vol[av_volfrac<=0.1], av_volfrac[av_volfrac<=0.1], linewidth=0, marker='.', markersize=10, color='#d95f0e')

pt.savefig('porosity_av_volfrac.pdf', format='pdf', bbox_inches = 'tight')
pt.savefig('porosity_av_volfrac.svg', bbox_inches = 'tight')

pt.clf()
pt.figure(figsize=(8, 2))
pt.xlabel(r'Cell volume [$\rm{\AA}^3$]')
pt.ylabel(r'Accessible volume [$cm^3/g$]')
pt.xlim([6000,9500])
pt.ylim([0,0.3])
pt.plot(uc_vol[av_cm3_g>0.1], av_cm3_g[av_cm3_g>0.1], linewidth=0, marker='.', markersize=10, color='#31a354')
pt.plot(uc_vol[av_cm3_g<=0.1], av_cm3_g[av_cm3_g<=0.1], linewidth=0, marker='.', markersize=10, color='#d95f0e')
pt.savefig('porosity_av_cm3_g.pdf', format='pdf', bbox_inches = 'tight')
pt.savefig('porosity_av_cm3_g.svg', bbox_inches = 'tight')
