#!/bin/sh
#
#PBS -N raspa
#PBS -q long
#PBS -l walltime=71:59:00
#PBS -l nodes=1:ppn=1
#

module load RASPA2


cd $PBS_O_WORKDIR

bash run_makegrid.sh ${fw} ${ff} ${nx} ${ny} ${nz} ${guest0}
