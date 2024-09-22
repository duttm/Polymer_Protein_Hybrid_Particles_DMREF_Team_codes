#!/bin/bash

#SBATCH --mail-user=mh1314@scarletmail.rutgers.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=nvt-lip
#SBATCH --partition RM-shared
#SBATCH -N 1 --ntasks-per-node=32
#SBATCH --time=03:00:00
#SBATCH --output=slurm.%N.%j.out
#SBATCH --export=ALL

GMX_APP=/ocean/projects/dmr170002p/hooten/dmref/gmx21_5.sif

$GMX_APP gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
$GMX_APP gmx mdrun -ntomp 32 -v -deffnm nvt
