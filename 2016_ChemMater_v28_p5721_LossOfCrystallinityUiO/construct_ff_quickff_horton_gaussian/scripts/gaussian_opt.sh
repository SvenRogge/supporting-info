#!/bin/sh
#
#PBS -N 12opt1
#PBS -q long 
#PBS -m e
#PBS -l walltime=72:00:00
#PBS -l nodes=1:ppn=24
#
## change to gent gaussian group

module load cluster
module load Gaussian/g09_D.01-intel-2015a-amd64-gpfs

export MKL_NUM_THREADS=1
export KMP_STACKSIZE=16m
export G09NOHOARD=1

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

g09 < gaussian.com > gaussian.log
formchk gaussian.chk gaussian.fchk && rm gaussian.chk
cp gaussian.* $ORIGDIR

echo "Contents of the directories after the calculation:"
echo "ORIGDIR:" $ORIGDIR
ls $ORIGDIR
echo "WORKDIR:" $WORKDIR
ls $WORKDIR

cd $ORIGDIR
rm -rf $WORKDIR


