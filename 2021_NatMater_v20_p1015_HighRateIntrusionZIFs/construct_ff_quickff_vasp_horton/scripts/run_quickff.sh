#!/bin/sh
#
#PBS -N ZIF8
#PBS -q long
#PBS -l walltime=72:00:00
#PBS -l nodes=1:ppn=all
#

module load QuickFF/2.2.0-intel-2015b-Python-2.7.10

ORIGDIR=$PBS_O_WORKDIR
WORKDIR=/local/$PBS_JOBID

mkdir -p $WORKDIR
cp $ORIGDIR/quickff_input.txt $WORKDIR
cp $ORIGDIR/quickff_sample.chk $WORKDIR

echo "Contents of the directories before the calculation:"
echo "ORIGDIR:" $ORIGDIR
ls
cd $WORKDIR
echo "WORKDIR:" $WORKDIR
ls

python qff.py -c quickff_input.txt quickff_sample.chk

echo "Contents of the directories after the calculation:"
echo "ORIGDIR:" $ORIGDIR
ls $ORIGDIR
echo "WORKDIR:" $WORKDIR
ls $WORKDIR

cp $WORKDIR/* $ORIGDIR
cd $ORIGDIR
rm -rf $WORKDIR


