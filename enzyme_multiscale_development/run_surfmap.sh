#!/bin/bash

# Using GRO file, SEQRES section of PDB file, and atom names and charges from 
# topology, run Surfmap on a CG config.

GROFILE_NOEXT=$1
ATOMSFILE=$2
CHARGESFILE=$3
SEQRESFILE=$4

SCRIPTDIR=/path/to/scripts

GROFILE=$GROFILE_NOEXT.gro
PDBFILE=$GROFILE_NOEXT.pdb
PQRFILE=$GROFILE_NOEXT.pqr

$SCRIPTDIR/gro2pdbcoords.py $GROFILE $SEQRESFILE > $PDBFILE

$SCRIPTDIR/gro2pqr.py $GROFILE $ATOMSFILE $CHARGESFILE

source /env_surfmap/bin/activate
surfmap -tomap electrostatics --png --keep --docker -pdb $PDBFILE -pqr $PQRFILE 
surfmap -tomap all --png --keep --docker -pdb $PDBFILE
