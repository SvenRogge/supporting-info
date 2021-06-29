import numpy as np
import h5py
from yaff import System, Cell
import matplotlib.pyplot as pt
from molmod.units import picosecond, angstrom


angles_closed = [61.83,171.30,188.7,298.17]
angles_open = [28.4,137.9,222.1,331.6]


def plot_hist(angles, xaxis, color, label):
    hist, bin_edges = np.histogram(angles, xaxis, density=True)
    pt.plot( 0.5*(bin_edges[1:] + bin_edges[:-1]), hist, color=color, label=label)
    pt.fill_between( 0.5*(bin_edges[1:] + bin_edges[:-1]), 0, hist, color=color, alpha=0.03)

loadings = ['%03d' %i for i in range(0,80,4)]

cmap = pt.cm.get_cmap('plasma_r')
colors = ['#377eb8','#4daf4a','#984ea3','#ff7f00']

pt.clf()
xaxis = np.arange(0,52,1)

for idx, loading in enumerate(loadings):
    angles1 = np.loadtxt('%s/dihedrals_1.txt' %loading)
    angles2 = np.loadtxt('%s/dihedrals_2.txt' %loading)
    angles3 = np.loadtxt('%s/dihedrals_3.txt' %loading)
    angles4 = np.loadtxt('%s/dihedrals_4.txt' %loading)
    plot_hist(np.hstack((angles2,angles3)), xaxis, color=cmap(1.*idx/20), label=loading)

pt.axvline(x=180-angles_closed[1], color='#e41a1c', linestyle='--', label='closed')
pt.axvline(x=180-angles_open[1], color='#a65628', linestyle='--', label='open')

pt.xlabel('Dihedral angle [deg]')
pt.ylabel('Histogram [1/deg]')
pt.xlim([0,50])
pt.ylim([0,0.06])

pt.legend(ncol=3)
pt.gcf().set_size_inches(9,4)
pt.savefig('dihedrals_hist_v3.svg',format='svg')
pt.savefig('dihedrals_hist_v3.pdf', bbox_inches='tight', file='pdf')
