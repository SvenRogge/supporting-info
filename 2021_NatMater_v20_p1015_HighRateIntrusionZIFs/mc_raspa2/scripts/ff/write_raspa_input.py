#!/usr/bin/env python

import os
import h5py as h5
import numpy as np
from glob import glob

from horton.io import IOData
from horton import Cell, log
log.set_level(log.silent)
from horton.periodic import periodic as pt
from horton.units import amu, angstrom, kjmol
from horton.periodic import periodic
from horton.constants import boltzmann
from yaff import System, log, ForceField
log.set_level(log.silent)

from cif import load_cif, dump_cif
from medff import prepare_medff, get_dispersion_coefficients

def write_mm3_mbis( fw = 'uio66-oien2014b'):
    workdir = os.path.join('..','ffpars','mm3-mbis', fw)
#    data = load_cif(os.path.join('..','cif','%s.cif'%fw))
    data, _ = prepare_medff(os.path.join('..','periodic',fw,'framework','PBE','mbis.h5'))
    system = System(data.numbers, data.pos, rvecs=data.cell.rvecs)
    system.detect_bonds()
#    print system.numbers[442]
#    print system.numbers[list(system.neighs1[442])]
    system.detect_ffatypes(ffatypes_mm3)
    # Read mbis charges for framework
    with h5.File(os.path.join('..','periodic',fw,'framework','PBE','mbis.h5'),'r') as fh5:
        charges = fh5['valence_charges'][:] + fh5['core_charges'][:]
    charges[:] -= np.mean(charges)
    iodata = IOData(numbers=system.numbers,coordinates=system.pos,cell=Cell(system.cell.rvecs))
    ffatypes = ["%s_%05d" % (system.ffatypes[system.ffatype_ids[iatom]],iatom) for iatom in xrange(system.natom)]
    dump_cif(os.path.join('..','ffpars','mm3-mbis',fw,'%s.cif'%fw), iodata, ffatypes=ffatypes)
    # Read mbis charges for guest molecule
    monomer, monomer_formula = 'methane','CH4'
    workdir = os.path.join('..','monomers',monomer,'B3LYP-aug-cc-pvtz')
    iodata = IOData.from_file(os.path.join(workdir,'gaussian.fchk'))
    with h5.File(os.path.join(workdir,'wpart_mbis.h5'),'r') as fh5:
        charges_monomer = fh5['valence_charges'][:] + fh5['core_charges'][:]
    charges_monomer[:] -= np.mean(charges)
    if monomer=='methane':
        ffatypes_monomer = np.array(['CCH4','HCH4','HCH4','HCH4','HCH4'])
        monomer_numbers = np.array([6,1,1,1,1])
    else: raise NotImplementedError
    ffas, q_monomer, z_monomer = [],[],[]
    for ffatype in np.unique(ffatypes_monomer):
        ffas.append(ffatype)
        q_monomer.append( np.mean(charges_monomer[ffatypes_monomer==ffatype]))
        z_monomer.append( monomer_numbers[ffatypes_monomer==ffatype][0] )
    # Write pseudo atoms file
    with open(os.path.join('..','ffpars','mm3-mbis',fw,'pseudo_atoms.def'),'w') as f:
        f.write("#number of pseudo atoms\n%d\n#type      print     as   chem  oxidation         mass       charge   polarization B-factor  radii  connectivity anisotropic anisotropic-type   tinker-type\n" % (system.natom+len(ffas)))
        for iatom in xrange(system.natom):
            ffa = system.get_ffatype(iatom) + ("_%05d" % iatom)
            f.write("%-12s yes %6s %6s %10d %12.8f %+32.20f         %6.3f   %6.3f %6.3f         %5d       %5d       %10s         %5d\n" % (ffa, pt[system.numbers[iatom]].symbol, pt[system.numbers[iatom]].symbol, 0, pt[system.numbers[iatom]].mass/amu, charges[iatom],0.0,1.0,1.0,0,0,"absolute",0) )
        for iatom in xrange(len(ffas)):
            f.write("%-12s yes %6s %6s %10d %12.8f %+32.20f         %6.3f   %6.3f %6.3f         %5d       %5d       %10s         %5d\n" % (ffas[iatom], pt[z_monomer[iatom]].symbol, pt[z_monomer[iatom]].symbol, 0, pt[z_monomer[iatom]].mass/amu, q_monomer[iatom],0.0,1.0,1.0,0,0,"absolute",0) )

def medff2raspa(system0, mbis0, upars, system1=None, mbis1=None, do_cross=True):
    '''
    Construct a MEDFF that RASPA can read based on MBIS parameters.
    See combination of RASPA documentation and MEDFF paper to understand where
    the formulas come from.
    '''
    if system1 is None:
        assert mbis1 is None
        system1, mbis1 = system0, mbis0
        do_cross = False
    else: do_cross = True and do_cross
    # Definition of interaction parameters
    U_olp = upars[0]+upars[1]
    U_s8 = upars[2]
    # Average the AIM parameters
    systems, mbiss, numberss = [], [], []
    for system, mbis in [(system0,mbis0),(system1,mbis1)]:
        newnumbers = []
        mbis_avg = {}
        for key in mbis.keys(): mbis_avg[key] = np.zeros(system.ffatypes.shape[0])
        for iffat, ffat in enumerate(system.ffatypes):
            mask = system.ffatype_ids==iffat
            number = system.numbers[mask][0]
            newnumbers.append(number)
            assert np.all(system.numbers[mask]==number)
            for key in mbis.keys():
                mbis_avg[key][iffat] = np.mean(mbis[key][mask])
        systems.append(system)
        mbiss.append(mbis_avg)
        numberss.append(np.array(newnumbers,dtype=int).copy())
    # Calculate resulting dispersion coefficients
    c6s, c8s, R_cross = get_dispersion_coefficients(
        np.concatenate((mbiss[0]['valence_widths'],mbiss[1]['valence_widths'])),
        np.concatenate((mbiss[0]['polar_scales'],mbiss[1]['polar_scales'])),
        np.concatenate((mbiss[0]['exp_r2'],mbiss[1]['exp_r2'])),
        np.concatenate((mbiss[0]['exp_r4'],mbiss[1]['exp_r4'])),
        np.concatenate((numberss[0],numberss[1])),)
    # Collect the atomic pair parameters
    nffa0 = len(systems[0].ffatypes)
    nffa1 = len(systems[1].ffatypes)
    data = []
    for iffa in xrange(nffa0):
        for jffa in xrange(nffa1):
            # No need to write same interaction twice
            if not do_cross and jffa>iffa: continue
            a,b = mbiss[0]['valence_widths'][iffa],mbiss[1]['valence_widths'][jffa]
            Nprod = mbiss[0]['valence_charges'][iffa]*mbiss[1]['valence_charges'][jffa]
            NaZb =  mbiss[0]['valence_charges'][iffa]*mbiss[1]['core_charges'][jffa]
            ZaNb =  mbiss[0]['core_charges'][iffa]*mbiss[1]['valence_charges'][jffa]
            # Electrostatic contribution
            p0 = -NaZb/boltzmann/angstrom
            p1 = -NaZb/a/2.0/boltzmann
            p2 = 0.0
            p3 = 0.0
            p4 = 0.0
            p5 = 1.0/a*angstrom
            p6 = -ZaNb/boltzmann/angstrom
            p7 = -ZaNb/b/2.0/boltzmann
            p8 = 1.0/b*angstrom
            # Use Taylor expansions
            if np.abs(a-b)<0.02:
                delta = b-a
                # Overlap contribution
                p0 += 0.0
                p1 += U_olp*Nprod*a**-3/384.0/np.pi*(6.0 - 9.0*delta/a)/boltzmann
                p2 += U_olp*Nprod*a**-4/384.0/np.pi*(6.0 - 9.0*delta/a)/boltzmann*angstrom
                p3 += U_olp*Nprod*a**-5/384.0/np.pi*(2.0 - 2.0*delta/a)/boltzmann*angstrom**2
                p4 += U_olp*Nprod*a**-6/384.0/np.pi*(1.0*delta/a)/boltzmann*angstrom**3
                # Electrostatic contribution
                p0 += -Nprod/boltzmann/angstrom
                p1 += -Nprod*(11.0/16.0 + 5.0*delta/32.0/a )/a/boltzmann
                p2 += -Nprod*(3.0/16.0 + 5.0*delta/32.0/a )/a/a/boltzmann*angstrom
                p3 += -Nprod*(1.0/48.0 + delta/16.0/a)/a/a/a/boltzmann*angstrom**2
                p4 += -Nprod*(delta/96.0/a)/a/a/a/a/boltzmann*angstrom**3
            else:
                # Overlap contribution
                p0 += U_olp*Nprod*4.0*a**2*b**2/(b**2-a**2)**3/boltzmann/angstrom/8.0/np.pi
                p1 += U_olp*Nprod*a/(b**2-a**2)**2/boltzmann/8.0/np.pi
                p6 += U_olp*Nprod*4.0*a**2*b**2/(a**2-b**2)**3/boltzmann/angstrom/8.0/np.pi
                p7 += U_olp*Nprod*b/(a**2-b**2)**2/boltzmann/8.0/np.pi
                # Electrostatic contribution
                p0 += -Nprod*a**4/(a**2-b**2)**2*(1.0-2*b**2/(a**2-b**2))/boltzmann/angstrom
                p1 += -Nprod*a**4/(a**2-b**2)**2/2.0/a/boltzmann
                p6 += -Nprod*b**4/(b**2-a**2)**2*(1.0-2*a**2/(b**2-a**2))/boltzmann/angstrom
                p7 += -Nprod*b**4/(b**2-a**2)**2/2.0/b/boltzmann

            # Dispersion parameters
            p9 = c6s[iffa,jffa+nffa0]/boltzmann/angstrom**6
            p10 = U_s8*c8s[iffa,jffa+nffa0]/boltzmann/angstrom**8
            p11 = R_cross[iffa,jffa+nffa0]*angstrom

            # Append to list
            data.append("%10s %10s %5s %+16.2f %+16.2f %12.2f %12.2f %12.2f %12.8f %16.2f %16.2f %12.8f %12.2f %12.2f %12.8f"%
                (systems[0].ffatypes[iffa],systems[1].ffatypes[jffa],"MEDFF",
            p0,p1,p2,p3,
            p4,p5,p6,p7,
            p8,p9,p10,p11))
    return data, systems, mbiss



def write_medff(fw='uio66-oien2014b',gg=None,suffix=''):
    '''
    Arguments:
        workdir:            path where files will be stores
        framework:          path to hdf5 file containing MBIS output and geometry for framework
        guest:              path to hdf5 file containing MBIS output and geometry for guest
        upars:              array containing interaction parameters [U_exch, U_ind, U_s8]
        upars_guestguest:   array containing interaction parameters [U_exch, U_ind, U_s8]
                            specifically for guest-guest interactions. If not provided, the
                            same upars as for guest-host interactions are used
        guest_ffatypes:     ATSELECT specification of atom types in the guest
                            If not provided, each atom is a specific atom type
    '''
    guests = ['methane','CO2']
    guest_ffatypes = [('CCH4','6&=4%1'),('HCH4','1&=1%6'),('C_co2','6&=2%8'),('O_co2','8&=1%(6&=2%8)')]
    # Load corresponding interaction parameters
    upars = {}
    fn_upars = '../medffpars/original.npy'
    upars['methane'] = np.load(fn_upars)
    fn_upars = '../medffpars/refit-ads-zif8-co2.npy'
    upars['CO2'] = np.load(fn_upars)
    print upars['CO2']
    print upars['methane']
    assert False
    # Construct directory to output
    ffname = 'medff%s'%suffix
    if gg is not None:
        ffname += "-ggua"
    workdir = os.path.join('..','ffpars',ffname,fw)
    if not os.path.isdir(workdir): os.makedirs(workdir)
    fn_upars = '../medffpars/original.npy'
    upars_guestguest = np.load(fn_upars)

    # Read the guest molecule
    system_guests, mbis_guests = [],[]
    for guest in guests:
        system_guest, mbis_guest = prepare_medff(os.path.join('..','monomers',guest,'B3LYP-aug-cc-pvtz','wpart_mbis.h5'), atypes=guest_ffatypes)
        system_guests.append(system_guest)
        mbis_guests.append(mbis_guest)
    # Read the framework
    system_framework, mbis_framework = prepare_medff(os.path.join('..','periodic',fw,'framework','PBE','mbis.h5'))
    # Write the framework geometry as a cif file
    iodata = IOData(numbers=system_framework.numbers, coordinates=system_framework.pos, cell=Cell(system_framework.cell.rvecs))
    dump_cif(os.path.join(workdir,'%s.cif'%fw), iodata, ffatypes=system_framework.ffatypes[system_framework.ffatype_ids])
    # Collect guest-guest interactions
    lines_ggs, lines_ghs, systemss, mbisss = [],[],[],[]
    for iguest, guest in enumerate(guests):
        for jguest in xrange(iguest,len(guests)):
            lines_gg,_,_ = medff2raspa(system_guests[iguest], mbis_guests[iguest], upars_guestguest, system1=system_guests[jguest], mbis1=mbis_guests[jguest],do_cross=iguest!=jguest)
            if gg is not None:
                if gg=='trappe' and guest=='methane' and guests[jguest]=='methane':
                    lines_gg = ["%10s %10s %10s %32.12f %32.12f" % ("CCH4","CCH4","LENNARD_JONES",147.90000000000000,3.730000000000000)]
                    lines_gg += ["%10s %10s %10s" % ("CCH4","HCH4","NONE")]
                    lines_gg += ["%10s %10s %10s" % ("HCH4","HCH4","NONE")]
                else: raise NotImplementedError
            lines_ggs += lines_gg
        # Collect guest-host interactions
        lines_gh, systems, mbiss = medff2raspa(system_guests[iguest], mbis_guests[iguest], upars[guest], system1=system_framework, mbis1=mbis_framework)
        if iguest==0:
            systemss.append(systems[1])
            mbisss.append(mbiss[1])
        lines_ghs += lines_gh
        systemss.append(systems[0])
        mbisss.append(mbiss[0])
#    lines_gh = medff2raspa(system_framework, mbis_framework, upars, system1=system_guest, mbis1=mbis_guest)
    # Turn off the host-host interactions
    unique_numbers = np.unique(system_framework.numbers)
    lines_hh = []
    for number0 in unique_numbers:
        for number1 in unique_numbers:
            if number0>number1: continue
            lines_hh.append("%9s_ %9s_ %5s" % (periodic[number0].symbol,periodic[number1].symbol,'NONE'))
    # Dump to file
    lines = lines_hh + lines_ghs + lines_ggs
    with open(os.path.join(workdir,'force_field.def'), 'w') as f:
        print >> f, "# rules to overwrite\n0\n# number of defined interactions\n%d\n# type type2  interaction"%(len(lines))
        for line in lines:
            print >> f, line
        print >> f, "# mixing rules to overwrite\n0"
    # Write a generic force_field_mixing_rules file, adapt to your needs!
    with open(os.path.join(workdir,'force_field_mixing_rules.def'),'w') as f:
        print >> f, "# general rule for shifted vs truncated\ntruncated\n# general rule tailcorrections\nyes\n# number of defined interactions\n0\n# type interaction\n# general mixing rule for Lennard-Jones\nLorentz-Berthelot"
    # Write the pseudo atoms file
    with open(os.path.join(workdir,'pseudo_atoms.def'),'w') as f:
        f.write("#number of pseudo atoms\n%d\n#type      print     as   chem  oxidation         mass       charge   polarization B-factor  radii  connectivity anisotropic anisotropic-type   tinker-type\n" % (np.sum([system.ffatypes.shape[0] for system in systemss])))
        for i in xrange(len(systemss)):
            system = systemss[i]
            mbis = mbisss[i]
            for iffa in xrange(system.ffatypes.shape[0]):
                iatom = np.where(system.ffatype_ids==iffa)[0][0]
                print >> f, "%-12s yes %6s %6s %10d %12.8f %+32.20f         %6.3f   %6.3f %6.3f         %5d       %5d       %10s         %5d" %\
                    (system.ffatypes[system.ffatype_ids[iatom]], periodic[system.numbers[iatom]].symbol, periodic[system.numbers[iatom]].symbol,\
                    0, periodic[system.numbers[iatom]].mass/amu, mbis['valence_charges'][iatom] + mbis['core_charges'][iatom],0.0,1.0,1.0,0,0,"absolute",0)

def write_saptff(fw):
    workdir = os.path.join('..','ffpars','saptff', fw)
    data = load_cif(os.path.join('..','cif','%s.cif'%fw))
    system = System(data['numbers'], data['coordinates'], rvecs=data['cell'].rvecs)
    system.detect_bonds()
    system.detect_ffatypes(ffatypes_saptff)
    iodata = IOData(numbers=system.numbers,coordinates=system.pos,cell=Cell(system.cell.rvecs))
    ffatypes = ["%s" % (system.ffatypes[system.ffatype_ids[iatom]]) for iatom in xrange(system.natom)]
    print os.path.join('..','ffpars','saptff',fw,'%s.cif'%fw)
    dump_cif(os.path.join('..','ffpars','saptff',fw,'%s.cif'%fw), iodata, ffatypes=ffatypes)


def write_saptff_general():
    solute, sfw, fw = read_saptff(write_raspa=True)
#    print solute, sfw, fw
    with open(os.path.join('..','ffpars','saptff','pseudo_atoms.def'),'w') as f:
        f.write("#number of pseudo atoms\n%d\n#type      print     as   chem  oxidation         mass       charge   polarization B-factor  radii  connectivity anisotropic anisotropic-type   tinker-type\n" % len(ffatypes_saptff))
        for iatom in xrange(len(ffatypes_saptff)):
            ffat = ffatypes_saptff[iatom][0]
            if '_' in ffat: symbol = ffat.split('_')[0]
            else: symbol = ffat[0]
            if ffat=='Zr': symbol='Zr'
            print ffat, symbol
            f.write("%-20s yes %6s %6s %10d %12.8f %+32.20f         %6.3f   %6.3f %6.3f         %5d       %5d       %10s         %5d\n" % (ffat, symbol, symbol, 0, periodic[symbol].mass/amu, 0.0,0.0,1.0,1.0,0,0,"absolute",0) )


ffatypes_dreiding =  [('Zn','30'),('O','8'),('C','6'),('N','7'),('Cl','17'),('H','1')]

def write_dreiding_mbis( fw = 'uio66-oien2014b', scheme='mbis'):
    workdir = os.path.join('..','ffpars','dreiding-%s-trappe'%scheme, fw)
    if not os.path.isdir(workdir): os.makedirs(workdir)
    data, _ = prepare_medff(os.path.join('..','periodic',fw,'framework','PBE','mbis.h5'))
    system = System(data.numbers, data.pos, rvecs=data.cell.rvecs)
    system.detect_bonds()
    system.detect_ffatypes(ffatypes_dreiding)
    # Read mbis charges for framework
    if scheme=='mbis':
        with h5.File(os.path.join('..','periodic',fw,'framework','PBE','mbis.h5'),'r') as fh5:
            charges = fh5['valence_charges'][:] + fh5['core_charges'][:]
    elif scheme in ['eqeq']:
        with h5.File(os.path.join('..','periodic',fw,'framework','eqeq','eqeq.h5'),'r') as fh5:
            charges = fh5['charges'][:]
    else: raise NotImplementedError
    charges[:] -= np.mean(charges)
    iodata = IOData(numbers=system.numbers,coordinates=system.pos,cell=Cell(system.cell.rvecs))
    ffatypes = ["%s_%05d" % (system.ffatypes[system.ffatype_ids[iatom]],iatom) for iatom in xrange(system.natom)]
    dump_cif(os.path.join('..','ffpars','dreiding-%s-trappe'%scheme,fw,'%s.cif'%fw), iodata, ffatypes=ffatypes)
    # Guest specifications
    ffatypes_monomer = np.array(['CH4_sp3','C_co2','O_co2','N_n2','N_com'])
    charges_monomer = np.array([0.0,0.70,-0.35,-0.482,0.964])
    monomer_numbers = np.array([6,6,8,7,7])
    monomer_masses = np.array([pt[6].mass + pt[1].mass*4,pt[6].mass,pt[8].mass,pt[7].mass,0.0])
    # Write pseudo atoms file
    with open(os.path.join('..','ffpars','dreiding-%s-trappe'%scheme,fw,'pseudo_atoms.def'),'w') as f:
        f.write("#number of pseudo atoms\n%d\n#type      print     as   chem  oxidation         mass       charge   polarization B-factor  radii  connectivity anisotropic anisotropic-type   tinker-type\n" % (system.natom+ffatypes_monomer.shape[0]))
        for iatom in xrange(system.natom):
            ffa = system.get_ffatype(iatom) + ("_%05d" % iatom)
            f.write("%-12s yes %6s %6s %10d %12.8f %+32.20f         %6.3f   %6.3f %6.3f         %5d       %5d       %10s         %5d\n" % (ffa, pt[system.numbers[iatom]].symbol, pt[system.numbers[iatom]].symbol, 0, pt[system.numbers[iatom]].mass/amu, charges[iatom],0.0,1.0,1.0,0,0,"absolute",0) )
        for iatom in xrange(ffatypes_monomer.shape[0]):
            f.write("%-12s yes %6s %6s %10d %12.8f %+32.20f         %6.3f   %6.3f %6.3f         %5d       %5d       %10s         %5d\n" % (ffatypes_monomer[iatom], pt[monomer_numbers[iatom]].symbol, pt[monomer_numbers[iatom]].symbol, 0, monomer_masses[iatom]/amu, charges_monomer[iatom],0.0,1.0,1.0,0,0,"absolute",0) )


def print_charges():
    data, _ = prepare_medff(os.path.join('..','periodic',fw,'framework','PBE','mbis.h5'))
    system = System(data.numbers, data.pos, rvecs=data.cell.rvecs)
    system.detect_bonds()
    ffatypes =  [('Zn','30'),('O','8'),('C1','6&=2%7'),('C3','6&=3%1'),('C2','6'),('N','7'),('Cl','17'),('H3','1&=1%(6&=3%1)'),('H1','1&=1%(6&=2%7)'),('H2','1')]
    system.detect_ffatypes(ffatypes)
    # Read mbis charges for framework
    with h5.File(os.path.join('..','periodic',fw,'framework','PBE','mbis.h5'),'r') as fh5:
        charges_mbis = fh5['valence_charges'][:] + fh5['core_charges'][:]
    charges_mbis[:] -= np.mean(charges_mbis)
    with h5.File(os.path.join('..','periodic',fw,'framework','eqeq','eqeq.h5'),'r') as fh5:
        charges_eqeq = fh5['charges'][:]
    charges_eqeq[:] -= np.mean(charges_eqeq)
    for iffa, ffa in enumerate(system.ffatypes):
        mask = np.where(system.ffatype_ids==iffa)[0]
        print "%5s %+8.2f %8.4f %+8.2f %8.4f" % (ffa, np.mean(charges_mbis[mask]), np.std(charges_mbis[mask]), np.mean(charges_eqeq[mask]), np.std(charges_eqeq[mask]))

def write_yaff(workdir, monomer='water'):
    system = System.from_file(os.path.join(workdir,'optcell.chk'))
    # Get the charges
    ff = ForceField.generate(system, os.path.join(workdir,'pars.txt'))

    iodata = IOData(numbers=system.numbers,coordinates=system.pos,cell=Cell(system.cell.rvecs))
    ffatypes = ["%s_%05d" % (system.ffatypes[system.ffatype_ids[iatom]],iatom) for iatom in xrange(system.natom)]
    fw = workdir.split(os.sep)[-1]		
    dump_cif(os.path.join(workdir,'%s.cif'%fw), iodata, ffatypes=ffatypes)


    charges = ff.system.charges
    if monomer=='water':
        ffas = ['O_TIP4P','H_TIP4P','M_TIP4P']
        z_monomer = np.array([8,1,0])
        q_monomer = np.array([0.0,0.5564,-2.0*0.5564])
        monomer_symbols = ['O','H','M']
        monomer_mass = [periodic[8].mass,periodic[1].mass,0.0]
    else: raise NotImplementedError


    # Write pseudo atoms file
    with open(os.path.join(workdir,'pseudo_atoms.def'),'w') as f:
        f.write("#number of pseudo atoms\n%d\n#type      print     as   chem  oxidation         mass       charge   polarization B-factor  radii  connectivity anisotropic anisotropic-type   tinker-type\n" % (system.natom+len(ffas)))
        for iatom in xrange(system.natom):
            ffa = system.get_ffatype(iatom) + ("_%05d" % iatom)
            f.write("%-12s yes %6s %6s %10d %12.8f %+32.20f         %6.3f   %6.3f %6.3f         %5d       %5d       %10s         %5d\n" % (ffa, pt[system.numbers[iatom]].symbol, pt[system.numbers[iatom]].symbol, 0, pt[system.numbers[iatom]].mass/amu, charges[iatom],0.0,1.0,1.0,0,0,"absolute",0) )
        for iatom in xrange(len(ffas)):
            f.write("%-12s yes %6s %6s %10d %12.8f %+32.20f         %6.3f   %6.3f %6.3f         %5d       %5d       %10s         %5d\n" % (ffas[iatom], monomer_symbols[iatom], monomer_symbols[iatom], 0, monomer_mass[iatom]/amu, q_monomer[iatom],0.0,1.0,1.0,0,0,"absolute",0) )

    print workdir

if __name__=='__main__':
    fns = sorted(glob(os.path.join('..','cif','zif8.cif')))
    for fn in fns:
        print fn
        fw = fn.split(os.sep)[-1].split('.')[0]
        write_dreiding_mbis(fw, scheme='mbis')

    ffname = 'qff-dreiding-mbis-tip4p'
    fws = [dn for dn in sorted(glob(os.path.join('..','ffpars',ffname,'*'))) if os.path.isdir(dn)]
    for fw in fws:
        write_yaff(fw)

