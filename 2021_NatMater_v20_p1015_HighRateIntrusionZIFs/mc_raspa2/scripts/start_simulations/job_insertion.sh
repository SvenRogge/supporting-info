#!/bin/sh
#
#PBS -N raspa
#PBS -q long
#PBS -l walltime=71:59:00
#PBS -l nodes=1:ppn=1
#

module load RASPA2/2.0.3-vtkrestart-intel-2018a


cd $PBS_O_WORKDIR

timeout 250000 bash run_insertion.sh ${fw} ${ff} ${nmol} ${nx} ${ny} ${nz} ${guest0}
if [ $? -eq 0 ]
then
    echo "Completed before timeout"
else
    echo "Rerun"
    qsub -N ${fw}_${nmol} -v fw=${fw},ff=${ff},nmol=${nmol},nx=${nx},ny=${ny},nz=${nz},guest0=${guest0} job_insertion.sh
fi
