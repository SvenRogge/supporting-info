#! /bin/bash

#PBS -N ElasticConstants_Cambridge
#PBS -q long
#PBS -l walltime=71:59:00
#PBS -l nodes=1:ppn=1
#

module load yaff/1.4.2-intel-2017b-Python-2.7.14

ORIGDIR=$PBS_O_WORKDIR
WORKDIR=/local/$PBS_JOBID

mkdir -p $WORKDIR

cp $ORIGDIR/* $WORKDIR/
cd $WORKDIR

python elastic_constants.py > output.log

cp $WORKDIR/* $ORIGDIR/
cd $ORIGDIR
rm -rf $WORKDIR


