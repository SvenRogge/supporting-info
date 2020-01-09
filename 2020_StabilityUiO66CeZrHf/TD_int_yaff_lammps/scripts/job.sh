#!/bin/bash
#
#PBS -N CeZr_2Zr4Hf1
#PBS -q long
#PBS -l walltime=72:00:00
#PBS -l nodes=1:ppn=4

source $VSC_DATA_VO/vsc40685_apps/activate.sh

module load yaff/1.4.2-intel-2017b-Python-2.7.14
module load LAMMPS/r12824-intel-2017b
module load mpi4py/2.0.0-intel-2017b-Python-2.7.14

ORIGDIR=$PBS_O_WORKDIR
WORKDIR=/local/$PBS_JOBID
SCRIPTDIR=$VSC_SCRATCH_VO/vsc40686/TD_int/UiO/ceriumMOFs_noupdate/2Zr4Hf1/scripts

mkdir -p $WORKDIR
cp $ORIGDIR/* $WORKDIR/.
cp $SCRIPTDIR/ymd.py $WORKDIR/.
cp $SCRIPTDIR/pars.txt $WORKDIR/.
cp $SCRIPTDIR/liblammps.py $WORKDIR/.
cp $SCRIPTDIR/mylammps.py $WORKDIR/.
cp $SCRIPTDIR/lammps_smoothei2.table $WORKDIR/.

cd $WORKDIR
mpirun -np 4 python ymd.py

cp $WORKDIR/* $ORIGDIR/.
cd $ORIGDIR
rm -rf $WORKDIR
rm $ORIGDIR/pars.txt
rm $ORIGDIR/ymd.py
rm $ORIGDIR/liblammps.py
rm $ORIGDIR/mylammps.py
rm $ORIGDIR/lammps_smoothei2.table

