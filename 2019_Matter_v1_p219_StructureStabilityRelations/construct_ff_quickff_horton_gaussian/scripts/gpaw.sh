#!/bin/sh
#
#PBS -N linker1_gpaw
#PBS -q long
#PBS -l walltime=72:00:00
#PBS -l nodes=1:ppn=20
#

. $VSC_DATA_VO/apps_cmm/activate.sh
module load gpaw/0.11.0.13004-intel-2015b-Python-2.7.10
module load gpaw-setups/0.9.11271-linux-x86_64

ORIGDIR=$PBS_O_WORKDIR
WORKDIR=/local/$PBS_JOBID

mkdir -p $WORKDIR
cp $ORIGDIR/* $WORKDIR

echo "Contents of the directories before the calculation:"
echo "ORIGDIR:" $ORIGDIR
ls
cd $WORKDIR
echo "WORKDIR:" $WORKDIR
ls

mpirun -np 20 gpaw-python ${HOME}/bin/gpaw-driver.py cluster.xyz --mode rg --spacing 0.15 --margin 4 --xc PBE --runtype energy 
 
echo "Contents of the directories after the calculation:"
echo "ORIGDIR:" $ORIGDIR
ls $ORIGDIR
echo "WORKDIR:" $WORKDIR
ls $WORKDIR

cp $WORKDIR/* $ORIGDIR
cd $ORIGDIR
rm -rf $WORKDIR


