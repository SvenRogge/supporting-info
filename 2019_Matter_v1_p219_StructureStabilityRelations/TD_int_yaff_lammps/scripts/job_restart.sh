#!/bin/bash
#
#PBS -N UiO_L6_TDint_r
#PBS -q long
#PBS -l walltime=72:00:00
#PBS -l nodes=1:ppn=4

source $VSC_DATA_VO/vsc40685_apps/activate.sh

module load yaff/1.1.3-intel-2017a-Python-2.7.13
module load LAMMPS/r12824-intel-2017a
module load mpi4py/2.0.0-intel-2017a-Python-2.7.13-timed-pingpong
module load vsc-mympirun/4.0.2

ORIGDIR=$PBS_O_WORKDIR
WORKDIR=/local/$PBS_JOBID

mkdir -p $WORKDIR
cp $ORIGDIR/* $WORKDIR/.

cd $WORKDIR
python add_cnt.py
mympirun python ymd_restart.py

cp $WORKDIR/* $ORIGDIR/.
cd $ORIGDIR
rm -rf $WORKDIR


