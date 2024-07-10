#!/bin/bash

# Using GRO file and SEQRES section of PDB file, run Surfmap on an AA config.


GROFILE_NOEXT=$1

ATOMSFILE=$2
CHARGESFILE=$3
SEQRES_FILE=$4

SCRIPTDIR=/path/to/scripts

GROFILE=$GROFILE_NOEXT.gro
PDBFILE=$GROFILE_NOEXT.pdb
PQRFILE=$GROFILE_NOEXT.pqr

$SCRIPTDIR/gro2pdbcoords.py $GROFILE $SEQRES_FILE > $PDBFILE


$SCRIPTDIR/gro2pqr.py $GROFILE $ATOMSFILE $CHARGESFILE

source /env_surfmap/bin/activate

export LD_LIBRARY_PATH=$HOME/apbs/APBS-3.4.1.Linux/lib:${LD_LIBRARY_PATH}
export PATH=$HOME/apbs/APBS-3.4.1.Linux/bin:${PATH}

surfmap -tomap electrostatics --png --keep --docker -pdb $PDBFILE -pqr $PQRFILE 
surfmap -tomap all --png --keep --docker -pdb $PDBFILE
