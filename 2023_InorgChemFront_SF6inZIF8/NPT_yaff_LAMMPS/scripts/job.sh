#!/bin/bash
#
#PBS -N ZIF8_6SF6
#PBS -l walltime=71:59:00
#PBS -l nodes=1:ppn=8

module load LAMMPS/7Aug2019-foss-2019b-Python-3.7.4-kokkos

ORIGDIR=$PBS_O_WORKDIR
WORKDIR=/local/$PBS_JOBID
mkdir -p $WORKDIR

cp $ORIGDIR/ymd.py $WORKDIR/.
cp $ORIGDIR/*.txt $WORKDIR/.
cp $ORIGDIR/*.chk $WORKDIR/.
cp $ORIGDIR/lammps.dat $WORKDIR/.
cp $ORIGDIR/lammps_smoothei2.table $WORKDIR/.

cd $WORKDIR
mpirun -np 8 python ymd.py

cp $WORKDIR/*.h5 $ORIGDIR/.

cd $ORIGDIR
rm -rf $WORKDIR
