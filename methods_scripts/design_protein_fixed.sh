#!/bin/bash

# Create Working Directory
WDIR=$SCRATCH/$JOB_NAME-$JOB_ID-$3
PROJECT_DIR=$WORK/test_allaa
DATABASE=$WORK/rosetta_src_2015.39.58186_bundle/main/database/

mkdir -p $WDIR
if [ ! -d $WDIR ]
then
  echo $WDIR not created
  exit
fi
cd $WDIR

# Copy Data and Config Files
cp $PROJECT_DIR/config/* .
cp $PROJECT_DIR/structures/$1.pdb .
# Put your Science related commands here
date
hostname

#./backrub -database $DATABASE -s $1.pdb -resfile NATAA.res -ex1 -ex2 -extrachi_cutoff 0 -backrub:mc_kt $2 -backrub:ntrials 10000 -nstruct 1 -out::suffix \_$2_$3 -backrub:initial_pack >& backrub.log.txt

date
./fixbb -database $DATABASE -s $1.pdb -resfile ALLAA.res -ex1 -ex2 -extrachi_cutoff 0 -nstruct 1 -overwrite -minimize_sidechains -linmem_ig 10 >& fixbb.log.txt

date

# Copy Results Back to Work Directory
RDIR=$PROJECT_DIR/results/$1/design-$1-$2-$JOB_ID-$3
mkdir -p $RDIR
cp * $RDIR/. 

# Cleanup 
rm -rf $WDIR
