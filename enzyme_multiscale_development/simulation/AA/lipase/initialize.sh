#!/bin/bash

GMX_APP=/home/mason/gmx21_5.sif

for i in {1..10}
do
# make a dir; cd to it
mkdir seed_$i; cd seed_$i
# copy in files up to solvation model including sim script
cp -a ../seed_0/* ./

#LOCAL - add ions and energy minimize
$GMX_APP gmx pdb2gmx -f 1oil.pdb -ff charmm36-feb2021 -water tip3p -o 1oil.gro
$GMX_APP gmx editconf -f 1oil.gro -c -d 1.5 -bt cubic -o 1oil_newbox.gro
$GMX_APP gmx solvate -cp 1oil_newbox.gro -cs spc216.gro -p topol.top -o 1oil_solv.gro

$GMX_APP gmx grompp -f ions.mdp -c 1oil_solv.gro -p topol.top -o ions.tpr
echo 15 | $GMX_APP gmx genion -s ions.tpr -p topol.top -pname NA -nname CL -neutral -o 1oil_solv_ions.gro

$GMX_APP gmx grompp -f minim.mdp -c 1oil_solv_ions.gro -p topol.top -o em.tpr
$GMX_APP gmx mdrun -v -deffnm em


#BRIDGES
# sbatch sim script
#sbatch simulation.bash
# grompp nvt
# mdrun nvt
# grompp md
# mdrun md

# return to dir up
cd ..
done
