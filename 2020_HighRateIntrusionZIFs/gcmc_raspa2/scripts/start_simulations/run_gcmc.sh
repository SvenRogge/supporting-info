#!/bin/bash

readarray pressures < pressures.dat
readarray temperatures < temperatures.dat

origdir=$(pwd)

framework=$1
forcefield=$2
temp=$3
press=$4
nx=$5
ny=$6
nz=$7
guest0=$8
dn=${framework}/${forcefield}/${guest0}/t${temp}-p${press}
# Find cut-off radius if present
if echo ${forcefield} | grep -w "rcut" > /dev/null
then
    rcut=${forcefield: -2}.0
else
    rcut=10.0
fi
printf "rcut = %s\n" $rcut
# Check whether electrostatics is present
chargemethod=ewald
coulombspacing=0.15
printf "Chargemethod = %s, coulombspacing = %s\n" $chargemethod $coulombspacing
if echo ${forcefield} | grep -w "ggua" > /dev/null
then
    omitaacoulomb=yes
else
    omitaacoulomb=no
fi
printf "Omitaacoulomb = %s\n" $omitaacoulomb


echo $dn
mkdir -p $dn
cd $dn
if [ -f simulation.input ]
then
  cp ${origdir}/restart.input restart.input 
  simulate -i restart.input >> raspa.log 2>> raspa.err
else
  cp ${origdir}/simulation.input simulation.input
  for fn in force_field.def force_field_mixing_rules.def ${guest0}.def ${guest1}.def pseudo_atoms.def ${framework}.cif
  do
      ln -sf ${origdir}/../ffpars/${forcefield}/${framework}/${fn} .
  done
  sed -i "s/__forcefield__/${forcefield}/g" simulation.input
  sed -i "s/__framework__/${framework}/g" simulation.input
  T=$(echo ${temperatures[$temp]} | tr '\n' ' ')
  sed -i "s/__temp__/${T}/g" simulation.input
  P=$(echo ${pressures[$press]} | tr '\n' ' ')
  sed -i "s/__press__/${P}/g" simulation.input
  sed -i "s/__nx__/${nx}/g" simulation.input
  sed -i "s/__ny__/${ny}/g" simulation.input
  sed -i "s/__nz__/${nz}/g" simulation.input
  sed -i "s/__guest0__/${guest0}/g" simulation.input
  sed -i "s/__rcut__/${rcut}/g" simulation.input
  sed -i "s/__chargemethod__/${chargemethod}/g" simulation.input
  sed -i "s/__coulombspacing__/${coulombspacing}/g" simulation.input
  sed -i "s/__omitaacoulomb__/${omitaacoulomb}/g" simulation.input
  
  if [ "$guest0" = "methane" ]
  then
      if echo ${forcefield} | grep -w "trappe" > /dev/null
      then
          ngrids=1
          guestatoms="CH4_sp3"
      else
          ngrids=2
          guestatoms="CCH4 HCH4"
      fi
      if [ "$forcefield" == "uff" ]
      then
          ngrids=2
          guestatoms="CCH4 HCH4"
      fi
  elif [ "$guest0" = "argon" ]
  then
      ngrids=1
      guestatoms=Ar
  elif [ "$guest0" = "N2" ]
  then
      ngrids=2
      guestatoms="N_n2 N_com"
  elif [ "$guest0" = "CO2" ]
  then
      ngrids=2
      guestatoms="C_co2 O_co2"
  elif [ "$guest0" = "water" ]
  then
      ngrids=3
      guestatoms="O_TIP4P H_TIP4P M_TIP4P"
  else
      exit 1
  fi
  echo $ngrids
  sed -i "s/__ngrids__/${ngrids}/g" simulation.input
  sed -i "s/__guestatoms__/${guestatoms}/g" simulation.input
  simulate -i simulation.input > raspa.log 2> raspa.err 
fi
cd $origdir
