#!/bin/bash
#
#PBS -N yaff_lammps_UiO66
#PBS -q long
#PBS -l walltime=72:00:00
#PBS -l nodes=1:ppn=4

source $VSC_DATA_VO/vsc40685_apps/activate.sh

module load yaff/1.0.develop.2.14-intel-2015b-Python-2.7.10-HDF5-1.8.15-patch1-serial
module load LAMMPS/r12824-intel-2015b
module load mpi4py/2.0.0-intel-2015b-Python-2.7.10
module load vsc-mympirun/3.4.2-intel-2015b-Python-2.7.10-vsc-base-2.4.2

ORIGDIR=$PBS_O_WORKDIR
WORKDIR=/local/$PBS_JOBID

mkdir -p $WORKDIR
cp $ORIGDIR/* $WORKDIR/.

cd $WORKDIR
mympirun python ymd.py

cp $WORKDIR/* $ORIGDIR/.
cd $ORIGDIR
rm -rf $WORKDIR


