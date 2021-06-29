#!/bin/bash
#
#PBS -N ZIF8_2-2-6_H2O_2-2-1_r1
#PBS -l walltime=72:00:00
#PBS -l nodes=1:ppn=18

source $VSC_DATA_VO/vsc40685_apps/activate.sh
module load LAMMPS/3Mar2020-intel-2019b-Python-3.7.4-kokkos

ORIGDIR=$PBS_O_WORKDIR
WORKDIR=/local/$PBS_JOBID
SCRIPTDIR1=$VSC_SCRATCH_VO/vsc40686/NPT_simulations/ZIF-8/waterIntrusion/water_gradient/scripts
SCRIPTDIR2=$VSC_SCRATCH_VO/vsc40686/NPT_simulations/ZIF-8/waterIntrusion/water_gradient/zif8_2-2-6_H2O_2-2-1/scripts
mkdir -p $WORKDIR

cp $ORIGDIR/restart_1.h5 $WORKDIR/.
cp $SCRIPTDIR2/init.chk $WORKDIR/.
cp $SCRIPTDIR1/pars.txt $WORKDIR/.
cp $SCRIPTDIR1/ghostatoms.py $WORKDIR/.
cp $SCRIPTDIR1/rcut_12.0.table $WORKDIR/.
cp $ORIGDIR/ymd_restart.py $WORKDIR/.

cd $WORKDIR
mpirun -np 18 python ymd_restart.py

cp $WORKDIR/* $ORIGDIR/.

cd $ORIGDIR
rm -rf $WORKDIR


