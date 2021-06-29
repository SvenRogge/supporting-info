#!/bin/bash
#
#PBS -N ZIF8_100K
#PBS -m ae
#PBS -l walltime=71:59:00
#PBS -l nodes=1:ppn=all

module purge
module load CP2K/5.1-intel-2018a
module load vsc-mympirun

cd ${PBS_O_WORKDIR}

ORIGDIR=$PBS_O_WORKDIR

WORKDIR=$VSC_SCRATCH_NODE/$PBS_JOBID

mkdir -p $WORKDIR
cd $WORKDIR

cp $ORIGDIR/* $WORKDIR/

(sleep 71h; touch EXIT_MD; sleep 20m; mkdir $ORIGDIR/temp; cp $WORKDIR/* $ORIGDIR/temp/ )&

#date
cd $WORKDIR
mympirun cp2k.popt -i input.inp -o output_1.out
#date

sleep 25m
kill $!

cp -r $WORKDIR/* $ORIGDIR
cd $ORIGDIR
rm -rf $WORKDIR
