#!/bin/bash
#
#PBS -N UiO66_12
#PBS -l walltime=71:59:00
#PBS -l nodes=1:ppn=18

module load LAMMPS/7Aug2019-foss-2019b-Python-3.7.4-kokkos

ORIGDIR=$PBS_O_WORKDIR
WORKDIR=/local/$PBS_JOBID
mkdir -p $WORKDIR

cp $ORIGDIR/init.chk $WORKDIR/.
cp $ORIGDIR/../../scripts/ymd.py $WORKDIR/.
cp $ORIGDIR/../input/*.txt $WORKDIR/
cp $ORIGDIR/../input/lammps_smoothei2.table $WORKDIR/.

cd $WORKDIR
mpirun -np 18 python ymd.py

cp $WORKDIR/*.h5 $ORIGDIR/.

cd $ORIGDIR
rm -rf $WORKDIR
