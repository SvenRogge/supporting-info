#!/usr/bin/env python

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import ticker
import matplotlib
matplotlib.rcParams.update({'font.size': 20, 'legend.fontsize':18})
from glob import glob
import subprocess

from horton import IOData, dump_cube, DenseLinalgFactory
from horton.grid import UniformGrid
from horton.constants import boltzmann
from horton.units import angstrom, kjmol
from horton import Cell
from yaff import System

from cif import load_cif

raspa_dir = '/home/steven/lib/RASPA2/share/raspa'

guestdirs={'CO2':'co2','methane':'ch4'}

def grep(fn, pattern):
    (out,err) = subprocess.Popen(['''grep "%s" %s | tail -n1'''%(pattern,fn)], stdout=subprocess.PIPE, shell=True).communicate()
    return out.split('\n')[:-1]

def write_cube_density(fw='ZIF8_opengate_111',ffname='qff-dreiding-mbis-tip4p',guest='water',nmol=20):
    basedir = os.path.join('..','raspa',fw,ffname,guest,'insertion_%d'%nmol)
    # Find number of particles per unit cell (density is already averaged over supercells)
    fn_out = sorted(glob(os.path.join(basedir,'Output','System_0','*.data')))
    assert len(fn_out)==1
    loading = float(grep(fn_out[0],'absolute adsorption')[0].split()[4][:-1])
    fn = os.path.join(basedir,'VTK','System_0','COMDensityProfile_%s.vtk'%guest)
    if not os.path.isfile(fn): return
    cycle = int(grep(fn_out[0],'cycle')[0].split()[2][:])
    fn_cube = os.path.join('..','densities','%s_%03d_%09d.cube'%(fw,nmol,cycle))
    if os.path.isfile(fn_cube):
        return
    print "%35s %8.2f %10d" % (fw, loading, cycle)
    fn_cif = os.path.join(basedir,'%s.cif'%fw)
    # Read the density from the VTK file
    def read_vtk(fn):
        with open(fn,'r') as f:
            for i in xrange(5):
                line = f.next()
            gridpoints = [int(w) for w in line.split()[1:]]
            density = np.zeros(gridpoints, dtype=float)
            for i in xrange(2):
                line = f.next()
            origin = np.array([float(w)*angstrom for w in line.split()[1:]])
            print origin
            for i in xrange(4):
                line = f.next()
            for i in xrange(density.shape[2]):
                for j in xrange(density.shape[1]):
                    for k in xrange(density.shape[0]):
                        density[k,j,i] = float(f.next())
            return origin, density
    origin, density = read_vtk(fn)
    if np.sum(density)==0.0:
        print "Density zero"
        return
    else:
        print "Sum density: %d" % (np.sum(density))
        print "Max density: %d" % (np.amax(density))
    density *= loading/np.sum(density)
    print origin, density.shape
    supercell = None
    npoints = np.array(density.shape)
    write_cube(fn_cif, fn_cube, density, origin=origin, reverse=True, label='Loading (cm$^3$(STP)/cm$^3$)', select_cross_sections=None, do_report=False, umax=None,direction='111',supercell=supercell,do_cross_sections=False)

def write_cube(fn_cif, fn, cubedata,umin=None,umax=None,unit=1.0,origin=np.zeros((3)),supercell=None,label='',reverse=False, do_cross_sections=True, select_cross_sections = None, do_report=False, direction='100', prefix='cubedata'):
    # Read the structure
    iodata = IOData.from_file(fn_cif)
    translation = np.zeros((3))
    translate_grid = np.array(np.array(cubedata.shape)*translation,dtype=int)
    translation = np.array(translate_grid,dtype=float)/np.array(cubedata.shape,dtype=float)
    cubedata = np.roll(np.roll(np.roll(cubedata, translate_grid[2],axis=2),translate_grid[1],axis=1),translate_grid[0],axis=0)
    # Make a supercell
    if supercell is not None:
        system = System(iodata.numbers, iodata.coordinates, rvecs=iodata.cell.rvecs)
        system = system.supercell(supercell[0],supercell[1],supercell[2])
        iodata = IOData(numbers=system.numbers, coordinates=system.pos, cell=Cell(system.cell.rvecs))

    # Put in central cell
    gvecs = iodata.cell.gvecs.ravel()
    center = np.zeros((3),dtype=long)
    for iatom in xrange(iodata.coordinates[:].shape[0]):
        cart = iodata.coordinates[iatom]
        center[0] = -np.ceil(gvecs[0]*cart[0] + gvecs[1]*cart[1] + gvecs[2]*cart[2] - 1.0)
        center[1] = -np.ceil(gvecs[3]*cart[0] + gvecs[4]*cart[1] + gvecs[5]*cart[2] - 1.0)
        center[2] = -np.ceil(gvecs[6]*cart[0] + gvecs[7]*cart[1] + gvecs[8]*cart[2] - 1.0)
        iodata.cell.add_rvec(iodata.coordinates[iatom], center)

    rvecs = np.zeros((3,3))
    for i in xrange(3):
        rvecs[i] = iodata.cell.rvecs[i]/cubedata.shape[i]#np.diag(br.delta)
    grid = UniformGrid(origin, rvecs, np.array(cubedata.shape), np.array([1,1,1]))
    cube = IOData(cube_data=cubedata/unit, coordinates=iodata.coordinates, numbers=iodata.numbers, grid=grid)
    dump_cube(fn,cube)
    extent = ((0.0,1.0,0.0,1.0))
    if umin is None: umin = np.amin(cubedata)
    if umax is None: umax = np.amax(cubedata)
    if do_cross_sections:
        if select_cross_sections is None: select_cross_sections = range(cubedata.shape[0])
        for xindex in select_cross_sections:
            if reverse: cmap = cm.seismic_r
            else: cmap = cm.seismic
            plt.clf()
            # Select data points along a given crystallographic direction
            if direction=='100':
                tdata = cubedata[xindex]
            elif direction=='111':
                n = cubedata.shape[0]
                assert np.all(np.array(cubedata.shape)==cubedata.shape[0])
                tdata = np.zeros((n,n))
                for ia in xrange(n):
                    for ib in xrange(n):
                        ic = (xindex - ia - ib)%n
                        tdata[ia,ib] = cubedata[ia,ib,ic]
            elif direction=='011':
                n = cubedata.shape[0]
                assert np.all(np.array(cubedata.shape)==cubedata.shape[0])
                tdata = np.zeros((n,n))
                for ia in xrange(n):
                    for ic in xrange(n):
                        ib = (xindex - ia)%n
                        tdata[ia,ic] = cubedata[ia,ia,ic]
            else: raise NotImplementedError

            im = plt.imshow(tdata/unit, interpolation='bilinear', origin='lower',
                        cmap=cmap, extent=extent, vmin=umin/unit, vmax=umax/unit)
            plt.xlabel("$x_b$")
            plt.ylabel("$x_c$")
            if not do_report:
                plt.title("$x_a=%4.2f$"%( (float(xindex)/cubedata.shape[0]) ), y=1.08)
                CBI = plt.colorbar(im, orientation='horizontal', shrink=0.8, label=label)
                tick_locator = ticker.MaxNLocator(nbins=6)
                CBI.locator = tick_locator
                CBI.update_ticks()
                CBI.solids.set_rasterized(True)
                CBI.solids.set_edgecolor("face")
            plt.locator_params(nbins=3)
            plt.gcf().set_size_inches(8,8)
            plt.tight_layout()
            if do_report:
                basename = '%s_report_%05d.svg'%(prefix,xindex)
            else:
                basename = '%s_%05d.png'%(prefix,xindex)
            plt.savefig(os.path.join(os.path.dirname(fn),basename))
            print "SAVED", os.path.join(os.path.dirname(fn),basename)
        if do_report:
            plt.clf()
            fig = plt.figure(figsize=(8, 1.5))
            ax1 = fig.add_axes([0.05, 0.80, 0.9, 0.15])
            norm = matplotlib.colors.Normalize(vmin=umin/unit, vmax=umax/unit)

#            # ColorbarBase derives from ScalarMappable and puts a colorbar
#            # in a specified axes, so it has everything needed for a
#            # standalone colorbar.  There are many more kwargs, but the
#            # following gives a basic continuous colorbar with ticks
#            # and labels.
            CBI = matplotlib.colorbar.ColorbarBase(ax1, cmap=cmap,
                                            norm=norm,
                                            orientation='horizontal', label=label)
            tick_locator = ticker.MaxNLocator(nbins=6)
            CBI.locator = tick_locator
            CBI.update_ticks()
            CBI.solids.set_rasterized(True)
            CBI.solids.set_edgecolor("face")
            plt.savefig(os.path.join(os.path.dirname(fn),'colorbar.svg'))
            plt.close('all')

def newgrid():
    fn_cif = os.path.join('..','cif','%s.cif'%'cok17')
    supercell=[2,2,2]
    iodata = load_cif(fn_cif, None)
    print iodata['cell'].rvecs/angstrom
    assert False

    system = System(iodata['numbers'], iodata['coordinates'], rvecs=iodata['cell'].rvecs)
    system = system.supercell(supercell[0],supercell[1],supercell[2])
    data = IOData(numbers=system.numbers, coordinates=system.pos, cell=Cell(system.cell.rvecs))

    nx = 100
    cubedata = np.zeros((nx,nx,nx))
    frac = [0.25,0.25,0.0]
    w = 5
    cubedata[int(frac[0]*nx)-w:int(frac[0]*nx)+w,int(frac[1]*nx)-w:int(frac[1]*nx)+w,0:int(frac[2]*nx)+w] = 1.0
    cubedata[int(frac[0]*nx)-w:int(frac[0]*nx)+w,int(frac[1]*nx)-w:int(frac[1]*nx)+w,int(frac[2]*nx)-w+nx:nx] = 1.0

    frac = [0.5,0.5,0.5]
    cubedata[int(frac[0]*nx)-w:int(frac[0]*nx)+w,int(frac[1]*nx)-w:int(frac[1]*nx)+w,int(frac[2]*nx)-w:int(frac[2]*nx)+w] = 2.0


    print "MAX", np.amax(cubedata)
    print int(frac[0]*nx)-w
    print int(frac[2]*nx)-w,int(frac[2]*nx)+w

    rvecs = data.cell.rvecs.copy()
    for i in xrange(3):
        rvecs[i] /= cubedata.shape[i]
    origin = np.zeros((3))
    grid = UniformGrid(origin, rvecs, np.array(cubedata.shape), np.array([1,1,1]))
    cube = IOData(cube_data=cubedata, coordinates=data.coordinates, numbers=data.numbers, grid=grid)
    dump_cube(fn_cif.replace('.cif','.cub'), cube)

if __name__=='__main__':
    for nmol in [20,40,60,80][:]:
        for suffix in ['opengate','closedgate'][:]:
            write_cube_density(fw='ZIF8_%s_111'%suffix,nmol=nmol)
            write_cube_density(fw='ZIF8_%s'%suffix,nmol=8*nmol)
