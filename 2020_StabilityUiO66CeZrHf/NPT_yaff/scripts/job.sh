#!/bin/bash
#
#PBS -N CeZr_NPT_xx
#PBS -q long
#PBS -l walltime=72:00:00
#PBS -l nodes=1:ppn=1

source $VSC_DATA_VO/vsc40685_apps/activate.sh
module load yaff/1.4.2-intel-2017b-Python-2.7.14

ORIGDIR=$PBS_O_WORKDIR
WORKDIR=/local/$PBS_JOBID
SCRIPTDIR=$VSC_SCRATCH_VO/vsc40686/NPT_simulations/UiO/ceriumMOFs_noupdate/xxx/scripts
mkdir -p $WORKDIR

cp $ORIGDIR/* $WORKDIR/.
cp $SCRIPTDIR/pars.txt $WORKDIR/.
cp $SCRIPTDIR/init.chk $WORKDIR

cd $WORKDIR
python ymd.py

cp $WORKDIR/* $ORIGDIR/.
rm $ORIGDIR/*.py
rm $ORIGDIR/pars.txt
rm $ORIGDIR/init.chk

cd $ORIGDIR
rm -rf $WORKDIR


