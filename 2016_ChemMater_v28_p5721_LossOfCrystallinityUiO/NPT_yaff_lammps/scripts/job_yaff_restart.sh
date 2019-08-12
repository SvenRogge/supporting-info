#!/bin/sh
#
#PBS -N UiO-66_def05
#PBS -q long
#PBS -l walltime=71:59:00
#PBS -l nodes=1:ppn=1

module load yaff/1.0.thesissven2.10.cmm-intel-2015b-Python-2.7.10
ORIGDIR=$PBS_O_WORKDIR
WORKDIR=/local/$PBS_JOBID

mkdir -p $WORKDIR
cp $ORIGDIR/restart_1.h5 $WORKDIR/.
cp $ORIGDIR/pars.txt $WORKDIR/.
cp $ORIGDIR/ymd_restart.py $WORKDIR/.

cd $WORKDIR
python ymd_restart.py

cp $WORKDIR/* $ORIGDIR/.
cd $ORIGDIR
rm -rf $WORKDIR

