#!/bin/sh

prefix='m'
suffix='MPa'
sign=''

#for i in 0MPa 100MPa 200MPa 300MPa 400MPa 500MPa 600MPa 700MPa 800MPa 900MPa 950MPa
for i in 0MPa 100MPa 900MPa
do
    cp ./ymd_restart.py ../${i}/
    cp ./add_cnt.py ../${i}
    cp ./job_restart.sh ../${i}/
    cp ./mylammps.py ../${i}
    cp ./liblammps.py ../${i}
    cp ./lammps_smoothei2.table ../${i}
    cp ./pars.txt ../${i}
done

#for i in 0MPa 100MPa 200MPa 300MPa 400MPa 500MPa 600MPa 700MPa 800MPa 900MPa 950MPa
for i in 0MPa 100MPa 900MPa
do
    cd ../${i}
    qsub job_restart.sh
done

