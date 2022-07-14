#!/bin/bash
#
#PBS -N job_id
#PBS -l walltime=71:59:00
#PBS -l nodes=1:ppn=1

module load LAMMPS/7Aug2019-foss-2019b-Python-3.7.4-kokkos

ORIGDIR=$PBS_O_WORKDIR
SCRIPTDIR=$PBS_O_WORKDIR/../scripts
WORKDIR=/local/$PBS_JOBID

mkdir -p $WORKDIR

cp $SCRIPTDIR/lammps_smoothei2.table $WORKDIR/.
cp $SCRIPTDIR/pars.txt $WORKDIR/.
cp $ORIGDIR/init.chk $WORKDIR/.
cp $ORIGDIR/ymd_us.py $WORKDIR/.

cd $WORKDIR
python ymd_us.py

cp $WORKDIR/*.h5 $ORIGDIR/.

cd $ORIGDIR
rm -rf $WORKDIR
