#!/usr/bin/env python

import os, sys
from glob import glob
import subprocess

def grep(fn, pattern):
    (out,err) = subprocess.Popen(['''grep "%s" %s'''%(pattern,fn)], stdout=subprocess.PIPE, shell=True).communicate()
    return out.split('\n')[:-1]
overwrite=True
if len(sys.argv)==1:
    pattern = '*/*/*'
else: pattern=sys.argv[1]
sims = sorted(glob('%s/t*-p*/Output/System_0/*.data'%pattern))
for sim in sims:
    comps=sim.split(os.sep)
    fw = comps[0]
    guest = comps[2]
    ff = comps[1]
    conds = comps[3]
    # Dump results to this file
    fn = sim.replace('.data','.summary')
    assert fn!=sim
    # If simulation output is not newer than summary, simply continue
    if not overwrite:
        if os.path.isfile(fn):
            if os.path.getmtime(sim)< os.path.getmtime(fn): continue
#    print sim
    # Read conversion factors
    out = grep(sim, "Conversion factor molecules")
    conversion = [float(line.split()[-2]) for line in out]
    conversions = [conversion[:len(conversion)/2],conversion[len(conversion)/2:]]
    # Read absolute loading
    out = grep(sim, "Average loading absolute \[molecules")
    index_cycle = 2
    loadings, errors = [], []
    skip = False
    print sim,out
    if len(out)==1:
        # Completed simulation
        loadings, errors = [], []
        for line in out:
            loading = float(line.split()[-4])
            if loading==0.0:
                error = 0.0
            else:
                error = float(line.split()[-2])/loading
            loadings.append(loading)
            errors.append(error)
    elif len(out)==0:
        # Simulation not yet completed
        out = grep(sim, "absolute adsorption")
        outs = [out[::2],out[1::2]]
        for out in outs:
            line = out[-1].split()
            if len(line)==8:
                # Equilibration not yet finished
                if len(out)>8:
                    loading = 0.0
                    for line in out[1:]:
                        loading += float(line.split()[2])
                    loading /= len(out)-1
                    index_cycle = 3
                    error = 0.0
                else:        
                    skip = True
                    continue
            elif len(line)==14:
                # Take average so far
                loading = float(line[4][:-1])
                error = 0.0
            else: raise NotImplementedError
            loadings.append(loading)
            errors.append(error)
    else: raise NotImplementedError
    if skip: continue
    # Read number of cycles completed
    out = grep(sim, "Current cycle")
    print sim, out[-1].split(), index_cycle
    cycles = int(out[-1].split()[index_cycle])
    # Read external pressure
    out = grep(sim, "External Pressure")
    pressure = float(out[-1].split()[2])
    #conds = sim.split(os.sep)
    fout = open(fn,'w')
    # Read molfractions
    molfrac = grep(sim, "MolFraction")
    molfracs = [float(line.split()[1]) for line in molfrac]
    print >> fout, "%30s %10s %30s %3s %3s %10d %20.5f" % (fw,guest, ff,conds.split('-')[0],conds.split('-')[1],cycles,pressure),
    for i in xrange(2):
        print >> fout, "%16.8f %16.8f" % (molfracs[i],loadings[i]),
        for conv in conversions[i]:
            print >> fout, "%16.8f" % (loadings[i]*conv),
        print >> fout, "%16.8f" % (errors[i]),
    print >> fout, "\n",
    fout.close()
