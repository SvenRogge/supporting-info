#!/bin/sh
#
#PBS -q long
#PBS -o output.file
#PBS -e error.file
#PBS -l walltime=30:05:00
#PBS -l nodes=1:ppn=4
#PBS -N mtd_test

module load yaff/1.0.develop.2.15-intel-2015b-Python-2.7.10-HDF5-1.8.15-patch1-serial
module load LAMMPS/r12824-intel-2015b
module load mpi4py/2.0.0-intel-2015b-Python-2.7.10
module load vsc-mympirun/3.4.2-intel-2015b-Python-2.7.10-vsc-base-2.4.2


ORIGDIR=$PBS_O_WORKDIR
WORKDIR=/local/$PBS_JOBID

mkdir -p $WORKDIR
rm -rf $ORIGDIR/data-lammps
mkdir $ORIGDIR/data-lammps
cp -r $ORIGDIR/* $WORKDIR/.

cd $WORKDIR


start=`date +%s`
python ff-lammps/gentable_pars.py
end=`date +%s`
runtime=$end-$start
echo "Runtime gentable: " $((runtime/3600)) "h" $((runtime%3600/60)) "m" $((runtime%60)) "s" > runtime.txt
start=`date +%s`
mympirun python m*.py 2000 1200 
end=`date +%s`
runtime=$end-$start
echo "Runtime md: " $((runtime/3600)) "h" $((runtime%3600/60)) "m" $((runtime%60)) "s" >> runtime.txt


cp -r $WORKDIR/* $ORIGDIR/.
cd $ORIGDIR
rm -rf $WORKDIR

