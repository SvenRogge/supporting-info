#!/usr/bin/env python

from glob import glob
import numpy as np
import math
import sys, os

from yaff import System, ForceField, ForcePart, PairPotMM3, Switch3, ForcePartPair, log, PairPotLJ
from molmod.units import angstrom, kjmol, pascal, kcalmol
mpa = 1e6*pascal

def corrections_mm3_switch3(epsilon, sigma, rc, w):
    cg,cg1,cg3,cg5 = epsilon,rc,2.0*sigma,w
    t1 = cg3 ** 2
    t2 = t1 * cg3
    t3 = 0.5e1 / 0.216e3 * t2
    t5 = cg1 / 9
    t10 = -cg1 + cg5
    t14 = t10 ** 2
    t17 = 0.1e1 / cg3
    t20 = math.exp(12 * t17 * cg5)
    t30 = cg1 ** 2
    t35 = math.exp(-12 * t17 * cg1)
    t37 = cg5 ** 2
    t39 = 0.1e1 / t37 / cg5
    t43 = cg1 * t10
    t44 = math.log(-t10)
    t47 = math.log(cg1)
    t54 = t1 ** 2
    t64 = cg * (0.6388888889e3 * ((-t3 + (0.7e1 / 0.36e2 * cg5 - t5) * t1 - 0.2e1 / 0.3e1 * (cg5 - cg1 / 4) * t10 * cg3 + cg5 * t14) * t20 + t3 + (cg5 / 12 + t5) * t1 + cg1 * (cg5 + cg1 / 3) * cg3 / 2 + t30 * cg5) * t35 * t2 * t39 - 0.225e1 * (2 * t43 * t44 - 2 * t43 * t47 + cg5 * (cg5 - 2 * cg1)) * t54 * t1 / cg1 / t10 * t39)
    # For a switching function, partial integration can be applied and the
    # following result holds
    return t64,-3.0*t64

def corrections_mm3(epsilon, sigma, rc):
    cg,cg1,cg3 = epsilon,rc,2.0*sigma
    t4 = math.exp(-12 / cg3 * cg1)
    t6 = cg1 ** 2
    t10 = cg3 ** 2
    t14 = t10 ** 2
    t21 = cg * (0.2129629630e3 * cg3 * t4 * (12 * cg1 * cg3 + t10 + 72 * t6) - 0.7500000000e0 * t14 * t10 / t6 / cg1)
    t4 = math.exp(-12 / cg3 * cg1)
    t5 = cg1 ** 2
    t6 = t5 * cg1
    t10 = cg3 ** 2
    t17 = t10 ** 2
    t23 = cg * (-0.6388888889e3 * t4 * (12 * cg1 * t10 + t10 * cg3 + 72 * t5 * cg3 + 288 * t6) + 0.450e1 * t17 * t10 / t6)
    return t21,t23

def corrections_lj_switch3(epsilon, sigma, rc,w):
    cg,cg1,cg3,cg5 = epsilon,rc,sigma,w
    t1 = cg3 ** 2
    t2 = t1 ** 2
    t3 = t2 * t1
    t5 = cg1 ** 2
    t6 = t5 * cg1
    t7 = t5 ** 2
    t8 = t7 * t6
    t9 = -cg1 + cg5
    t10 = t9 ** 2
    t12 = t10 ** 2
    t13 = t12 * t10 * t9
    t14 = t8 * t13
    t15 = math.log(-t9)
    t18 = math.log(cg1)
    t24 = cg5 ** 2
    t25 = t24 ** 2
    t26 = t25 * t24
    t31 = t25 * cg5
    t37 = t7 ** 2
    t43 = t24 * cg5
    t60 = t26 * t7 * t5 - t26 * t3 / 84 - 6 * t31 * t8 + t31 * cg1 * t3 / 18 + 15 * t25 * t37 - t25 * t5 * t3 / 9 - 20 * t43 * t37 * cg1 + t43 * t6 * t3 / 9 + 15 * t24 * t37 * t5 - t24 * t7 * t3 / 18 - 6 * t37 * t6 * cg5 + t37 * t7
    t71 = -4 * t3 * cg * (2 * t14 * t15 - 2 * t14 * t18 + cg5 * (cg5 - 2 * cg1) * t60) / t8 / t43 / t13
    return t71, -3.0*t71

def corrections_lj(epsilon, sigma, rc):
    t1 = cg3 ** 2
    t2 = t1 ** 2
    t3 = t2 * t1
    t5 = cg1 ** 2
    t6 = t5 ** 2
    t10 = t6 ** 2
    t16 = -0.4e1 / 0.9e1 * t3 * cg * (3 * t6 * t5 - t3) / t10 / cg1
    t1 = cg3 ** 2
    t2 = t1 ** 2
    t3 = t2 * t1
    t5 = cg1 ** 2
    t6 = t5 ** 2
    t11 = t6 ** 2
    t17 = 0.8e1 / 0.3e1 * t3 * cg * (3 * t6 * t5 - 2 * t3) / t11 / cg1
    return t16, t17


class ForcePartTailCorrection(ForcePart):
    def __init__(self, system, pair_pot):
        super(ForcePartTailCorrection, self).__init__('pair_tailcorr', system)
        self.ecorr = 0.0
        self.wcorr = 0.0
        self.system = system
        assert isinstance(pair_pot, PairPotMM3) or isinstance(pair_pot, PairPotLJ)
        results = {}
        for iffa in xrange(system.ffatypes.shape[0]):
            for jffa in xrange(system.ffatypes.shape[0]):
                i = np.where(system.ffatype_ids==iffa)[0][0]
                j = np.where(system.ffatype_ids==jffa)[0][0]
                sigma = 0.5*(pair_pot.sigmas[i] + pair_pot.sigmas[j])
                epsilon = np.sqrt(pair_pot.epsilons[i]*pair_pot.epsilons[j])
                tr = pair_pot.get_truncation()

                if isinstance(pair_pot, PairPotLJ):
                    if isinstance(tr,Switch3):
                        e,w = corrections_lj_switch3(epsilon,sigma,pair_pot.rcut,pair_pot.get_truncation().width)
                    elif tr is None:
                        e,w = corrections_lj(epsilon,sigma,pair_pot.rcut)
                    else:
                        raise NotImplementedError
                elif isinstance(pair_pot, PairPotMM3):
                    if isinstance(tr,Switch3):
                        e,w = corrections_mm3_switch3(epsilon,sigma,pair_pot.rcut,pair_pot.get_truncation().width)
                    elif tr is None:
                        e,w = corrections_mm3(epsilon,sigma,pair_pot.rcut)
                    else:
                        raise NotImplementedError
                else: raise NotImplementedError
                imask = system.ffatype_ids==iffa
                jmask = system.ffatype_ids==jffa
                fac = np.sum(imask)*np.sum(jmask)
                self.ecorr += 2.0*np.pi*e*fac
                self.wcorr += 2.0*np.pi*w*fac

    def _internal_compute(self, gpos, vtens):
        if self.system.cell.nvec != 3: return 0.0
        if vtens is not None:
            w = self.wcorr/self.system.cell.volume/3.0
            vtens[0,0] = w
            vtens[1,1] = w
            vtens[2,2] = w
        return self.ecorr/self.system.cell.volume
