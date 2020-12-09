#! /usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as pt
import numpy as np

def process_xrd(loading):
    fns = ['xrd_traj%d.txt' %i for i in range(1,21)]
    data_0 = np.loadtxt('%s/%s' %(loading,fns[0]))
    x_axis = data_0[:,0]
    intensity_total = data_0[:,1]

    for fn in fns[1:]:
        data = np.loadtxt('%s/%s' %(loading,fn))
        intensity_total += data[:,1]

    return x_axis, intensity_total/len(fns)


loadings = ['%03d' %i for i in range(0,80,4)]

cmap = pt.cm.get_cmap('plasma_r')

pt.clf()

for idx, loading in enumerate(loadings):
    avg_x, avg_y = process_xrd(loading)
    pt.plot(avg_x, avg_y+40*idx, color=cmap(1.*idx/20), label=loading)
    
pt.xlabel(r"2$\theta$ (deg)")
pt.ylabel("Intensity")
pt.xlim([0,50])
lgd = pt.legend(bbox_to_anchor=(1.05, 1.0), ncol=2)

pt.gca().xaxis.grid(True)


pt.tick_params(
    axis='y',          # changes apply to the y-axis
    which='both',      # both major and minor ticks are affected
    left='off',      # ticks along the left edge are off
    right='off',         # ticks along the right edge are off
    labelleft='off') # labels along the left edge are off

bbox_extra_artists=(lgd,)
pt.savefig("pxrd_all_v2.pdf", bbox_extra_artists=(lgd,), bbox_inches='tight', file='pdf')
pt.savefig("pxrd_all_v2.svg", bbox_extra_artists=(lgd,), bbox_inches='tight', file='svg')
