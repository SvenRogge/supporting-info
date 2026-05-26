#!/usr/bin/env python

import os,sys,copy,numpy as np
from yaff import *
from molmod import MolecularGraph

def get_indices(bonds, numbers):
    graph = MolecularGraph(bonds, numbers)
    indices = graph.independent_vertices
    return indices

def translate(P,v):
    trans_matrix= np.array( [[1,0,0,v[0]],[0,1,0,v[1]],[0,0,1,v[2]],[0,0,0,1]])
    todo = np.ones(4)
    todo[0] = P[0]
    todo[1] = P[1]
    todo[2] = P[2]
    return np.dot(trans_matrix,todo)

def get_cv(ff):
    ls = get_indices(ff.system.bonds, ff.system.numbers)
    l1 = ls[1]
    l2 = ls[2]
    layers = [l1,l2]
    return [CVCOMProjection(ff.system,groups=layers,index=0),CVCOMProjection(ff.system,groups=layers,index=1)]

def adapt_structure(ff,cv):
    cv_object = get_cv(ff)
    # Adapt stucture for cv
    pos = copy.copy(ff.system.pos)
    ls = get_indices(ff.system.bonds, ff.system.numbers)
    for idx in ls[1]:
        pos[idx] = translate(pos[idx],[*cv,0])[0:3]
    ff.update_pos(pos)
    cell_symmetrize(ff)

    try:
        assert np.allclose(np.array(cv)[:2],np.array([cv_object[0].compute(), cv_object[1].compute()]))
    except AssertionError:
        print(cv, np.array([cv_object[0].compute(), cv_object[1].compute()]))
        sys.exit(1)
