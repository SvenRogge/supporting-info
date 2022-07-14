#!/bin/sh
#
#PBS -N qff_SF6
#PBS -l walltime=0:59:00
#PBS -l nodes=1:ppn=1
#

ORIGDIR=$PBS_O_WORKDIR
WORKDIR=/local/$PBS_JOBID

mkdir -p $WORKDIR
cp $ORIGDIR/SF6_quickff_sample_2.chk $WORKDIR
cp $ORIGDIR/quickff_input.txt $WORKDIR
cp $ORIGDIR/pars_ei.txt $WORKDIR
cp $ORIGDIR/pars_lj.txt $WORKDIR

echo "Contents of the directories before the calculation:"
echo "ORIGDIR:" $ORIGDIR
ls
cd $WORKDIR
echo "WORKDIR:" $WORKDIR
ls

source $VSC_DATA_VO/apps_cmm/activate.sh
module load yaff/1.6.0-intel-2020a-Python-3.8.2
module load QuickFF/2.2.7-intel-2020a-Python-3.8.2

# Run covalent QuickFF ff derivation
qff.py -c quickff_input.txt SF6_quickff_sample_2.chk 

echo "Contents of the directories after the calculation:"
echo "ORIGDIR:" $ORIGDIR
ls $ORIGDIR
echo "WORKDIR:" $WORKDIR
ls $WORKDIR

cp $WORKDIR/* $ORIGDIR
cd $ORIGDIR
rm -rf $WORKDIR
