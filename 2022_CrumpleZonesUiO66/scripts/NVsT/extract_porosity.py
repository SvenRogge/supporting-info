#! /usr/bin/env python

from os import path

with open('porosity.csv', 'w') as g:
    g.write('#unit_cell_volume,density,AV_A3,AV_volumefraction,AV_cm3g\n')

for i in range(6000, 9825, 25):

    if path.isfile('%s/struc_av.vol' %i):
        with open('%s/struc_av.vol' %i, 'r') as f:
            line = f.readline()
            line_split = line.split()
            
        uc_vol = line_split[3]
        density = line_split[5]
        av_a3 = line_split[7]
        av_volfrac = line_split[9]
        av_cm3_g = line_split[11]
        
        with open('porosity.csv', 'a') as g:
            g.write('%s,%s,%s,%s,%s\n' %(uc_vol,density,av_a3,av_volfrac,av_cm3_g))
