#!/bin/bash
#PBS -N /user/home/gent/vsc409/vsc40944/doc/periodic/ZIF8/V0/rerun/high/rerun/rerun/rerun
#PBS -m ae
#PBS -q long
#PBS -l walltime=71:59:00
#PBS -l nodes=1:ppn=all

module purge
export VSC_INSTITUTE_CLUSTER=golett
ulimit -s unlimited

module load cluster

module load jobs
module load scripts
module load VASP/5.4.1-intel-2015b-mt-vaspsol-20150914
module list

cd /user/home/gent/vsc409/vsc40944/doc/periodic/ZIF8/V0/rerun/high/rerun/rerun/rerun
free -m
date
touch CHGCAR WAVECAR

rm /user/home/gent/vsc409/vsc40944/doc/periodic/ZIF8/V0/rerun/high/rerun/rerun/rerun/tempout
mympirun --output /user/home/gent/vsc409/vsc40944/doc/periodic/ZIF8/V0/rerun/high/rerun/rerun/rerun/tempout vasp
free -m
date

exit 0
