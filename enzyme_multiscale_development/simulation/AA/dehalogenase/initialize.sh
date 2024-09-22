#!/bin/bash

GMX_APP=/home/mason/gmx21_5.sif
ROOT_DIR=$PWD
for i in {1..10}
do
cd $ROOT_DIR
# make a dir; cd to it
mkdir seed_$i; cd seed_$i
# copy in files up to solvation model including sim script
cp -a ../seed_0/* ./

#LOCAL - add ions and energy minimize
$GMX_APP gmx pdb2gmx -f 3rk4-CL.pdb -ff charmm36-feb2021 -water tip3p -o 3rk4.gro
$GMX_APP gmx editconf -f 3rk4.gro -c -d 1.5 -bt cubic -o 3rk4_newbox.gro
$GMX_APP gmx solvate -cp 3rk4_newbox.gro -cs spc216.gro -p topol.top -o 3rk4_solv.gro

$GMX_APP gmx grompp -f ions.mdp -c 3rk4_solv.gro -p topol.top -o ions.tpr
echo 15 | $GMX_APP gmx genion -s ions.tpr -p topol.top -pname NA -nname CL -neutral -o 3rk4_solv_ions.gro

$GMX_APP gmx grompp -f minim.mdp -c 3rk4_solv_ions.gro -p topol.top -o em.tpr
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
