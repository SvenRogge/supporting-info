#!/bin/bash

for dn in $(ls -d ZIF8_*gate/*/*/t*p*)
do
    grep "absolute adsorption" ${dn}/Output/System_0/*.data | grep avg | awk -F ' ' '{print $3" " $5}' | sed 's/.$//' > ${dn}/gcmc.out
    ncycle=$(grep "cycle" ${dn}/Output/System_0/*.data | tail -n1 | cut -d " " -f 3)
    adsorption=$(grep "absolute adsorption" ${dn}/Output/System_0/*.data | tail -n1)
    printf "%75s %12d Current: %10s Average: %10s\n" $dn $ncycle ${adsorption:22:10} ${adsorption:37:10}
done
