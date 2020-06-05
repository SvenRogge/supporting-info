#!/bin/sh
#
#PBS -N 12_gpaw
#PBS -q long
#PBS -l walltime=72:00:00
#PBS -l vmem=400GB
#PBS -l nodes=1:ppn=12
#

module load gpaw/0.11.0.13004-intel-2015b-Python-2.7.10
module load gpaw-setups/0.9.11271-linux-x86_64
module load scripts

ORIGDIR=$PBS_O_WORKDIR
WORKDIR=/local/$PBS_JOBID

mkdir -p $WORKDIR
cp $ORIGDIR/cluster.xyz $WORKDIR

echo "Contents of the directories before the calculation:"
echo "ORIGDIR:" $ORIGDIR
ls
cd $WORKDIR
echo "WORKDIR:" $WORKDIR
ls

mympirun gpaw-python ${HOME}/bin/gpaw-driver.py cluster.xyz --mode rg --spacing 0.15 --xc PBE --margin 4 --runtype sp

echo "Contents of the directories after the calculation:"
echo "ORIGDIR:" $ORIGDIR
ls $ORIGDIR
echo "WORKDIR:" $WORKDIR
ls $WORKDIR

cp $WORKDIR/* $ORIGDIR
cd $ORIGDIR
rm -rf $WORKDIR


