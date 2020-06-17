#!/bin/bash

origdir=$(pwd)

framework=$1
forcefield=$2
nx=$3
ny=$4
nz=$5
guest=$6
dn=${framework}/${forcefield}/${guest}/makegrid
# Find cut-off radius if present
if echo ${forcefield} | grep -w "rcut" > /dev/null
then
    rcut=${forcefield: -2}.0
else
    rcut=12.0
fi
printf "rcut = %s\n" $rcut
# Check whether electrostatics is present
if ( echo ${forcefield} | grep -w "trappe" > /dev/null || [ "${forcefield}" = "uff" ] || [ "${forcefield}" = "saptff" ] ) && [ "$guest" = "methane" ]
then
    chargemethod=none
    coulombspacing=10.0
else
    chargemethod=ewald
    coulombspacing=0.15
fi
printf "Chargemethod = %s, coulombspacing = %s\n" $chargemethod $coulombspacing
# Check if adsorbate adsorbate guest guest interactions need to be included
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
cp ${origdir}/makegrid.input simulation.input
for fn in force_field.def force_field_mixing_rules.def ${guest}.def pseudo_atoms.def ${framework}.cif
do
    ln -sf ${origdir}/../ffpars/${forcefield}/${framework}/${fn} .
done
sed -i "s/__forcefield__/${forcefield}/g" simulation.input
sed -i "s/__framework__/${framework}/g" simulation.input
sed -i "s/__nx__/${nx}/g" simulation.input
sed -i "s/__ny__/${ny}/g" simulation.input
sed -i "s/__nz__/${nz}/g" simulation.input
sed -i "s/__guest__/${guest}/g" simulation.input
sed -i "s/__rcut__/${rcut}/g" simulation.input
sed -i "s/__chargemethod__/${chargemethod}/g" simulation.input
sed -i "s/__coulombspacing__/${coulombspacing}/g" simulation.input
sed -i "s/__omitaacoulomb__/${omitaacoulomb}/g" simulation.input
if [ "$guest" = "methane" ]
then
#    if echo "uff-trappe dreiding-uff-trappe" | grep -w $forcefield > /dev/null
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
elif [ "$guest" = "argon" ]
then
    ngrids=1
    guestatoms=Ar
elif [ "$guest" = "N2" ]
then
    ngrids=2
    guestatoms="N_n2 N_com"
elif [ "$guest" = "CO2" ]
then
    ngrids=2
    guestatoms="C_co2 O_co2"
elif [ "$guest" = "water" ]
then
    ngrids=3
    guestatoms="O_TIP4P H_TIP4P M_TIP4P"
else
    exit 1
fi
sed -i "s/__ngrids__/${ngrids}/g" simulation.input
sed -i "s/__guestatoms__/${guestatoms}/g" simulation.input
simulate 2> raspa.log
cd $origdir

