#!/bin/sh
#
#PBS -N raspa
#PBS -l walltime=71:59:00
#PBS -l nodes=1:ppn=1
#

module load RASPA2


cd $PBS_O_WORKDIR
#fw=zif8-jelle-prim-md-40-08000
#ff=qff-dreiding-mbis-tip4p
#t=0
#p=0
#nx=2
#ny=2
#nz=2
#guest0=water

timeout 250000 bash run_gcmc.sh ${fw} ${ff} ${t} ${p} ${nx} ${ny} ${nz} ${guest0}
if [ $? -eq 0 ]
then
    echo "Completed before timeout"
else
    echo "Rerun"
    qsub -N ${fw}_${p}_${guest0} -v fw=${fw},ff=${ff},t=${t},p=${p},nx=${nx},ny=${ny},nz=${nz},guest0=${guest0} job_gcmc.sh
fi
