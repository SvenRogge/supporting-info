#!/bin/sh
#
#PBS -N def_type1
#PBS -q long
#PBS -l walltime=71:59:00
#PBS -l nodes=1:ppn=1

module load /apps/gent/SL6/sandybridge/modules/all/molmod/1.1-intel-2015b-Python-2.7.10
module load yaff/1.0.thesissven2.10.cmm-intel-2015b-Python-2.7.10

ORIGDIR=$PBS_O_WORKDIR
WORKDIR=/local/$PBS_JOBID

mkdir -p $WORKDIR
cp $ORIGDIR/init.chk $WORKDIR/.
cp $ORIGDIR/pars.txt $WORKDIR/.
cp $ORIGDIR/ymd.py $WORKDIR/.

cd $WORKDIR
unset LC_CTYPE
python ymd.py

cp $WORKDIR/* $ORIGDIR/.
cd $ORIGDIR
rm -rf $WORKDIR

