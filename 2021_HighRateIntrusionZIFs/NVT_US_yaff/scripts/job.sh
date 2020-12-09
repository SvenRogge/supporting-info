#!/bin/bash
#
#PBS -N job_id
#PBS -l walltime=71:59:00
#PBS -l nodes=1:ppn=1

source $VSC_DATA_VO/vsc40685_apps/activate.sh
source $VSC_DATA_VO/vsc41948_apps/activate.sh

module purge
module load yaff/1.6.0-intel-2019b-Python-3.7.4
module load LAMMPS/patch_20Nov2019-intel-2019b

ORIGDIR=$PBS_O_WORKDIR
SCRIPTDIR=$PBS_O_WORKDIR/../../inputfiles
WORKDIR=/local/$PBS_JOBID

mkdir -p $WORKDIR
cp $SCRIPTDIR/rcut_*.table $WORKDIR
cp $SCRIPTDIR/ghostatoms.py $WORKDIR
cp $SCRIPTDIR/pars*.txt $WORKDIR
cp $SCRIPTDIR/liblammps.py $WORKDIR
cp $ORIGDIR/* $WORKDIR/.

cd $WORKDIR
python ymd_ghost.py > yaff.log

cp $WORKDIR/* $ORIGDIR/.
rm $ORIGDIR/pars*.txt
rm $ORIGDIR/*.table
rm $ORIGDIR/lammps*.dat
rm $ORIGDIR/ghostatoms.py
rm $ORIGDIR/liblammps.py
cd $ORIGDIR
rm -rf $WORKDIR
