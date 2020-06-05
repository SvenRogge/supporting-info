#! /bin/bash

#PBS -N 12-Quickff
#PBS -q long
#PBS -l walltime=72:00:00
#PBS -l nodes=1:ppn=2
#

module load VSC-tools/1.7.1-ictce-4.1.13-scoop
module load quickff/1.1.1-qref-ictce-4.1.13-Python-2.7.3

ORIGDIR=$PBS_O_WORKDIR
WORKDIR=/local/$PBS_JOBID

mkdir -p $WORKDIR
cp $ORIGDIR/sample.chk $WORKDIR
cd $WORKDIR

qffest=$(which qff-est.py)

myscoop -h 2 $qffest --scoop --ei-model=Zero --vdw-model=Zero --ic-ids=all --suffix=_mbis --fn-traj=traj.pp sample.chk > quickff.log

cp $WORKDIR/* $ORIGDIR
cd $ORIGDIR
rm -rf $WORKDIR

