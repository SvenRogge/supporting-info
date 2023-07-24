#!/bin/sh
#
#PBS -N SF6
#PBS -l walltime=03:00:00
#PBS -l nodes=1:ppn=1
#
## change to gent gaussian group

ml purge
module load Gaussian/g16_A.03-intel-2019a

ORIGDIR=$PBS_O_WORKDIR
WORKDIR=/local/$PBS_JOBID

mkdir -p $WORKDIR
cp $ORIGDIR/gaussian.com $WORKDIR

echo "Contents of the directories before the calculation:"
echo "ORIGDIR:" $ORIGDIR
ls
cd $WORKDIR
echo "WORKDIR:" $WORKDIR
ls

g16 < gaussian.com > gaussian.log
formchk gaussian.chk gaussian.fchk && rm gaussian.chk
cp gaussian.* $ORIGDIR

echo "Contents of the directories after the calculation:"
echo "ORIGDIR:" $ORIGDIR
ls $ORIGDIR
echo "WORKDIR:" $WORKDIR
ls $WORKDIR

cd $ORIGDIR
rm -rf $WORKDIR

