#!/bin/sh
#
#PBS -N Fe4N
#PBS -q long
#PBS -l walltime=71:59:00
#PBS -l nodes=1:ppn=1
#

ORIGDIR=$PBS_O_WORKDIR
WORKDIR=/local/$PBS_JOBID

mkdir -p $WORKDIR
cp $ORIGDIR/sp.gpw $WORKDIR

echo "Contents of the directories before the calculation:"
echo "ORIGDIR:" $ORIGDIR
ls
cd $WORKDIR
echo "WORKDIR:" $WORKDIR
ls

# Set base directory, only delcatty supported for now
BASEDIR=/user/data/gent/gvo000/gvo00003/vsc40685_apps/delcatty
# Make the python packages visible
export PYTHONPATH=${BASEDIR}/lib/python/:$PYTHONPATH
# Make the executables visible
export PATH=${BASEDIR}/bin:$PATH
# Make the libraries visible
export LD_LIBRARY_PATH=${BASEDIR}/include:$LD_LIBRARY_PATH

# Load some modules for horton
module load libxc/2.2.2-intel-2015b
module load Libint/2.0.3-intel-2015b
module load h5py/2.5.0-intel-2015b-Python-2.7.10-serial
module load gpaw/0.11.0.13004-intel-2015b-Python-2.7.10
module load gpaw-setups/0.9.11271-linux-x86_64

# Convert CHGAR to denspart format (specify pseudo numbers!)
denspart-convert-gpaw sp.gpw denspart.h5
# Run MBIS partitioning
denspart-mbis denspart.h5 mbis.h5

echo "Contents of the directories after the calculation:"
echo "ORIGDIR:" $ORIGDIR
ls $ORIGDIR
echo "WORKDIR:" $WORKDIR
ls $WORKDIR

cp $WORKDIR/* $ORIGDIR
cd $ORIGDIR
rm -rf $WORKDIR

