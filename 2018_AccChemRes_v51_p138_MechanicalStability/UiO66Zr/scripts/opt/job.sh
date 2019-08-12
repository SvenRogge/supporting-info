#!/bin/sh
#
#PBS -N UiO-66_opt
#PBS -q long
#PBS -l walltime=71:59:00
#PBS -l nodes=1:ppn=1

module load yaff/1.0.thesissven2.11.cmm-intel-2016a-Python-2.7.11
ORIGDIR=$PBS_O_WORKDIR
WORKDIR=/local/$PBS_JOBID

mkdir -p $WORKDIR
cp $ORIGDIR/opt.chk $WORKDIR/.
cp $ORIGDIR/pars.txt $WORKDIR/.
cp $ORIGDIR/StrainStress.py $WORKDIR/.

cd $WORKDIR
python StrainStress.py

cp $WORKDIR/* $ORIGDIR/.
cd $ORIGDIR
rm -rf $WORKDIR

