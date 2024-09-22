#!/bin/bash

# Run this script in the production root directory
# This script initialize2.sh should follow initialize1.sh!

#REFERENCES
#http://md.chem.rug.nl/index.php/2021-martini-online-workshop/tutorials/564-2-proteins-basic-and-martinize-2#martini3-proteins

GMX_APP=/home/mason/gmx21_5.sif
AAREF_DIR=/home/mason/dutt_lab/dmref/lipase/aa/production
#GMX_APP=/ocean/projects/dmr170002p/hooten/dmref/gmx21_5.sif
#AAREF_DIR=/ocean/projects/dmr170002p/hooten/dmref/lipase/aa/production
ROOT_DIR=$PWD # Run in Production root directory

for i in {1..10}
do
cd $ROOT_DIR
cd seed_$i

# Convert to gro and fix box size same as AA ref
BOXDIMS=`tail lipaa_seed${i}_10ns.gro -n 1`
$GMX_APP gmx editconf -f lipase_cg.pdb -box $BOXDIMS -o lipase_cg.gro

# em in vacuum. pilots showed this was fine, tutorial agrees
$GMX_APP gmx grompp -f minim.mdp -c lipase_cg.gro -p topol_cg.top -o em.tpr
$GMX_APP gmx mdrun -pin on -v -deffnm em

## solvate. larger radius is from tutorial
$GMX_APP gmx solvate -cp em.gro -cs water.gro -radius 0.21 -o lipase_cg_sol.gro

########################### MANUAL STEP ###########################
# open topol_cg.top
# Add the water molecules to [ molecules ]
# save as topol_cg.top
########################################################################

done
