#!/bin/sh
#
#PBS -N mbis_horton
#PBS -l walltime=11:59:00
#PBS -l nodes=1:ppn=1
#

module load horton/2.1.1-intel-2020a-Python-2.7.18

ORIGDIR=$PBS_O_WORKDIR
WORKDIR=/local/$PBS_JOBID
mkdir -p $WORKDIR
cp $ORIGDIR/gaussian.fchk $WORKDIR

echo "Contents of the directories before the calculation:"
echo "ORIGDIR:" $ORIGDIR
ls

cd $WORKDIR
echo "WORKDIR:" $WORKDIR
ls

horton-wpart.py --grid=veryfine gaussian.fchk gaussian_wpart.h5:mbis mbis > horton_mbis.log

cp $WORKDIR/* $ORIGDIR

echo "Contents of the directories after the calculation:"
echo "ORIGDIR:" $ORIGDIR
ls $ORIGDIR
echo "WORKDIR:" $WORKDIR
ls $WORKDIR

cd $ORIGDIR
rm -rf $WORKDIR


