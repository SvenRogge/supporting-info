#!/bin/bash
#
#PBS -N ZIF8_grad_r
#PBS -q long
#PBS -l walltime=72:00:00
#PBS -l nodes=1:ppn=1

source $VSC_DATA_VO/vsc40685_apps/activate.sh
module load yaff/1.4.2-intel-2017b-Python-3.6.3

ORIGDIR=$PBS_O_WORKDIR
WORKDIR=/local/$PBS_JOBID
SCRIPTDIR=$VSC_SCRATCH_VO/vsc40686/NPT_simulations/ZIF-8/waterIntrusion/gradient/scripts
mkdir -p $WORKDIR

cp $ORIGDIR/init.chk $WORKDIR/.
cp $ORIGDIR/restart_19.h5 $WORKDIR/.
cp $SCRIPTDIR/pars.txt $WORKDIR/.
cp $SCRIPTDIR/ghostatoms.py $WORKDIR/.
cp $SCRIPTDIR/tailcorr.py $WORKDIR/.
cp $SCRIPTDIR/ymd_restart.py $WORKDIR/.

cd $WORKDIR
python ymd_restart.py

cp $WORKDIR/* $ORIGDIR/.
rm $ORIGDIR/*.py
rm $ORIGDIR/pars.txt

cd $ORIGDIR
rm -rf $WORKDIR


