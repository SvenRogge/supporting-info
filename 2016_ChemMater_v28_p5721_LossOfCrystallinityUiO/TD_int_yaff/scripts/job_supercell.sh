#!/bin/sh
#
#PBS -N UiO-66
#PBS -q long
#PBS -l walltime=71:59:00
#PBS -l nodes=1:ppn=1

module load /apps/gent/SL6/sandybridge/modules/all/molmod/1.1-intel-2015b-Python-2.7.10
module load yaff/1.0.develop.2.14-intel-2015b-Python-2.7.10-HDF5-1.8.15-patch1-serial
ORIGDIR=$PBS_O_WORKDIR
WORKDIR=/local/$PBS_JOBID

mkdir -p $WORKDIR
cp $ORIGDIR/init.chk $WORKDIR/.
cp $ORIGDIR/pars.txt $WORKDIR/.
cp $ORIGDIR/ymd_sc.py $WORKDIR/.

cd $WORKDIR
python ymd_sc.py

cp $WORKDIR/* $ORIGDIR/.
cd $ORIGDIR
rm -rf $WORKDIR

