#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import os
from glob import glob

from horton.units import angstrom
from horton import IOData

def plot_densities(sims, fig_fn, labels=None, nslabs=5, plotslabs=None, rhomin=1e-8, rhomax=1e-3):
    '''
    Plot 2D cross sections along z-axis of a density from a cube file.

    sims:
        list of tuples (fw,nmol) where it is assumed that the cube file is
        located in ../densities/fw_nmol.cube
    fig_gn:
        where to save the plot
    labels:
        list of len(sims), naming each simulations
    nslabs:
        number of slabs in the z-direction
    plotslabs:
        which slabs to plot, defaults to all slabs
    rhomin, rhomax:
        thresholds for minimal and maximal density
    '''
    densities = []
    for fn,nmol in sims:
        #fn = os.path.join('..','densities','%s_%03d.cube'%(fw,nmol))
        if not os.path.isfile(fn):
            print "WARNING, could not find file %s" % fn
            return
        iodata = IOData.from_file(fn)
        # Only tested for cubic cells, probably works for other too but the division
        # into slabs will be less intuitive
        L = iodata.cell.rvecs[0,0]
        assert np.allclose(iodata.cell.rvecs,L*np.eye(3))
        # Density is in number of molecules per unit cell, you might consider
        # converting this to molecules per volume
        density = iodata.cube_data
        densities.append(density)
    # Divide into slabs in z-direction
    if plotslabs is None: plotslabs = range(nslabs)
    assert density.shape[2]%nslabs==0
    slabsize = density.shape[2]/nslabs
    for symmetrize in [True, False]:
        plt.clf()
        # Color map settings
        cmap = plt.get_cmap('jet')
        cmap.set_under(color='black')
        cmap.set_over(color='grey')
        for isim, density in enumerate(densities):
            if symmetrize:
                # Symmetrize over equivalent pores
                for s in density.shape: assert s%2==0
                shift = tuple([s/2 for s in density.shape])
                density = 0.5*(density+np.roll(np.roll(np.roll(density,shift[0],axis=0),shift[1],axis=1),shift[2],axis=2))
            if labels is not None:
                plt.subplot(len(plotslabs),len(sims),isim+1)
                plt.title(labels[isim])
            for islab,slab in enumerate(plotslabs):
                data = np.sum(density[:,:,slab*slabsize:(slab+1)*slabsize],axis=2)
                plt.subplot(len(plotslabs),len(sims),len(sims)*islab+1+isim)
                print slab*slabsize,(slab+1)*slabsize
                plt.imshow(data,interpolation='spline36',cmap=cmap,vmin=rhomin,vmax=rhomax)
                xlabelnames = [0.0,0.5,1.0]
                xlabels = [data.shape[0]*x for x in xlabelnames]
                plt.xticks(xlabels,xlabelnames)
                plt.yticks(xlabels,xlabelnames)
        plt.gcf().set_size_inches(5*len(sims),4*len(plotslabs))
        if symmetrize:
            plt.savefig('../plots/%s_symm.png'%(fig_fn),bbox_inches="tight")
        else:
            plt.savefig('../plots/%s.png'%(fig_fn),bbox_inches="tight")

if __name__=='__main__':
    suffixes = ['opengate','closedgate']
    ucs = [('111','_111'),('222','')]
    for nmol in [4,8,20,40,60,80][:]:
    #for nmol in [4,8][:]:
        for uc, ucname in ucs:
            ucsize = np.prod([int(d) for d in uc])
            print ucsize
            sims = []
            for suffix in suffixes:
                print '../densities/ZIF8_%s%s_%03d_*'%(suffix,ucname,nmol)
                fns = sorted(glob('../densities/ZIF8_%s%s_%03d_*'%(suffix,ucname,nmol*ucsize)))
                if len(fns)==0: continue
                # Study convergence
                conv = [(fn,nmol*ucsize) for fn in fns]
                labels = ["$N=%5.1e$" % (float(fn.split('_')[-1].split('.')[0])) for fn in fns]
                plot_densities( conv, fig_fn='convergence_%s_%s_%03d'%(suffix,uc,nmol), labels=labels)
                sims.append( (fns[-1],nmol) )
            if len(sims)==0: continue
            plot_densities(sims, fig_fn = 'compare_densities_%03d_%s'%(nmol,uc), labels=suffixes)
            #assert False
