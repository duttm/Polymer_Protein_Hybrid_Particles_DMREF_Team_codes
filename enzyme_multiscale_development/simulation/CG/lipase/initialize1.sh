#!/bin/bash

# Run this script in the production root directory
# 

#REFERENCES
#http://md.chem.rug.nl/index.php/2021-martini-online-workshop/tutorials/564-2-proteins-basic-and-martinize-2#martini3-proteins

GMX_APP=/home/mason/gmx21_5.sif
AAREF_DIR=/home/mason/dutt_lab/dmref/lipase/aa/production
#GMX_APP=/ocean/projects/dmr170002p/hooten/dmref/gmx21_5.sif
#AAREF_DIR=/ocean/projects/dmr170002p/hooten/dmref/lipase/aa/production
ROOT_DIR=$PWD

for i in {1..10}
do
cd $ROOT_DIR
mkdir seed_$i; cd seed_$i

# each AA ref sim provides the initial structure for one CG prod sim 
cp -a $AAREF_DIR/seed_$i/lipaa_seed${i}_10ns.gro ./

# Martinize to obtain a protein-only PDB file and topology
martinize2 -f lipaa_seed${i}_10ns.gro -x lipase_cg-proteinonly.pdb -ff martini3001 -from universal -p backbone -ignh -ignore SOL -ignore NA -ignore CAL -elastic
find . -type f -name "molecule_*.itp" -exec mv {} {}.backup \;

# template_seed_0 contains common files
#martini_v3.0.0_solvents_v1.itp
#martini_v3.0.0.itp
#martini_v3.0.0_ions_v1.itp
#md.mdp
#minim.mdp
#molecule_0.itp
#nvt1.mdp
#nvt2.mdp
#nvt3.mdp
#simulation.bash
#topol_cg.top
#water.gro
cp -a ../template_seed_0/* ./

########################### MANUAL STEP ###########################
# Open lipase_cg-proteinonly.pdb
# STEP ONE
# Put the ION coords back in the PDB file with Martini types
#  e.g.
#   AA gro file lines like:
#  321CAL    CAL 4635   3.335   4.953   4.933
#30803NA      NA96079   0.481   2.292   7.804
#30804NA      NA96080   2.995   0.174   5.812
#   becomes CG pdb lines like:
#ATOM    702 CA  ION    321      33.35   49.53   49.33   1.00  0.00
#ATOM    703 NA  ION    322      04.81   22.92   78.04   1.00  0.00
#ATOM    704 NA  ION    322      29.95   01.74   58.12   1.00  0.00
# Save as lipase_cg.pdb
####
#NOTE : ion_gro_to_pdb.py helps automate this process
####
# STEP TWO
# Look for and remove erroneous TER lines within the protein...
########################################################################

########################### IGNORE STEP ###########################
# This file is identical for all seeds prior to CG solvation.

# Open topol_cg-proteinonly.top
# Add IONs to topol_cg-proteinonly.top
# Update include lines for correct force field and topologies
#include "martini_v3.0.0.itp"
#include "martini_v3.0.0_ions_v1.itp"
#include "martini_v3.0.0_solvents_v1.itp"
#include "molecule_0.itp"
# Update ion counts under [ molecules ]
#molecule_0    1
#CA 1
#NA 2
# Save as topol_cg.top
########################################################################

########################### IGNORE STEP ###########################
# Martinize does not reliably create the correct bonds when using a gro input.
# This is probably due to the reference structure crossing PBC boundaries,
# so it should be easy to fix by using trjconv -center -pbc mol (tested, works).
# The file lipase_cg.itp copied from seed_0 is a correct protein topology
# created using the martinize command above and modified as outlined below.
# Since this part of the config is the same in all seeds, 
# this file can be prepared in advance and reused.
# Its reference structure and therefore its elastic network are based on 
# AA ref seed 1 at 10ns.
# MODIFICATION
## Open molecule_0.itp
## In the lipase model, Martinize does apply the same protonation state to
# HIS204 and HIS286 that CHARMM does. Correct the SC3 bead on those residues to match the martini3/aminoacids.ff/HIH definition. Also change the charge from 0 to 1. 
## Looks like this in a corrected file
##448 TQ2p 204 HIS SC3 448  1.0
##629 TQ2p 286 HIS SC3 629  1.0
## Save as lipase_cg.itp.
########################################################################

done
