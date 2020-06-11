#!/bin/sh

prefix='m'
suffix='MPa'
sign=''

#for i in 0MPa 100MPa 200MPa 300MPa 400MPa 500MPa 600MPa 700MPa 800MPa 900MPa 1000MPa 1100MPa 1200MPa
for i in 950MPa 1050MPa
do
    if [[ $i == ${prefix}* ]]; then 
        sign='-'
        press=${i#$prefix}
    else
        sign=''
        press=${i}
    fi

    press=${press%$suffix}

    mkdir ../${i}
    cp ./ymd.py ../${i}/
    sed -i "s/press=1e6\*pascal/press=${sign}${press}e6*pascal/" ../${i}/ymd.py
    cp ./job.sh ../${i}/
    cp ./init.chk ../${i}/
    cp ./mylammps.py ../${i}
    cp ./liblammps.py ../${i}
    cp ./lammps_smoothei2.table ../${i}
    cp ./pars.txt ../${i}
done

#for i in 0MPa 100MPa 200MPa 300MPa 400MPa 500MPa 600MPa 700MPa 800MPa 900MPa 1000MPa 1100MPa 1200MPa
for i in 950MPa 1050MPa
do
    cd ../${i}
    qsub job.sh
done

