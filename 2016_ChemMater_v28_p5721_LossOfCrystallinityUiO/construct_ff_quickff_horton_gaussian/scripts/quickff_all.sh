#! /bin/bash

#PBS -N 12-Quickff
#PBS -q long
#PBS -l walltime=72:00:00
#PBS -l nodes=1:ppn=8
#

module load VSC-tools/1.7.1-ictce-4.1.13-scoop
module load quickff/1.1.1-qref-ictce-4.1.13-Python-2.7.3

ORIGDIR=$PBS_O_WORKDIR
WORKDIR=/local/$PBS_JOBID

mkdir -p $WORKDIR
cp $ORIGDIR/sample.chk $WORKDIR
cp $ORIGDIR/gaussian_wpart.h5 $WORKDIR
cd $WORKDIR

qffest=$(which qff-est.py)

myscoop -h 8 $qffest --scoop --ei-model=HarmGauss --ei-scales=1.0,1.0,1.0 --ei-path=mbis --vdw-model=Zero --ic-ids=all --suffix=_mbis --fn-traj=traj.pp sample.chk gaussian_wpart.h5 > quickff.log

cp $WORKDIR/* $ORIGDIR
cd $ORIGDIR
rm -rf $WORKDIR

