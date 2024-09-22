#!/bin/bash

# Run this script in the production root directory
# This script initialize3.sh should follow initialize2.sh!

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

## em in solvent
$GMX_APP gmx grompp -f minim.mdp -c lipase_cg_sol.gro -p topol_cg_sol.top -o em_solv.tpr
$GMX_APP gmx mdrun -v -deffnm em_solv

done

echo "initial setup complete! equilibration and md in simulation script."
