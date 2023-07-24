#!/bin/bash
#
#PBS -N UiO66_12_av
#PBS -l walltime=71:59:00
#PBS -l nodes=1:ppn=8

source $VSC_DATA_VO/vsc40685_apps/activate.sh
module load LAMMPS/7Aug2019-foss-2019b-Python-3.7.4-kokkos

ORIGDIR=$PBS_O_WORKDIR
WORKDIR=/local/$PBS_JOBID
mkdir -p $WORKDIR

cp $ORIGDIR/ymd.py $WORKDIR/.
cp $ORIGDIR/../input/*.txt $WORKDIR/
cp $ORIGDIR/../input/*.chk $WORKDIR/

cd $WORKDIR
mpirun -np 8 python ymd.py

cp $WORKDIR/*.h5 $ORIGDIR/.

cd $ORIGDIR
rm -rf $WORKDIR
