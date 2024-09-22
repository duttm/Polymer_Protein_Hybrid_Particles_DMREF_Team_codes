#!/bin/bash

#SBATCH --mail-user=mh1314@scarletmail.rutgers.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=lipcg
#SBATCH --partition=RM-shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=8:00:00
#SBATCH --output=slurm.%N.%j.out      # STDOUT output file
#SBATCH --export=ALL               # Export you current env to the job env

GMX_APP=/ocean/projects/dmr170002p/hooten/dmref/gmx21_5.sif

# restrained NVT equlibration -- give the water time to increase its entropy / get into the protein
$GMX_APP gmx grompp -f nvt1.mdp -c em_solv.gro -r em_solv.gro -p topol_cg_sol.top -o nvt1.tpr
$GMX_APP gmx mdrun -ntomp 16 -v -deffnm nvt1

# unrestrained NVT equlibration -- allow the system to relax into the box ~lento~
$GMX_APP gmx grompp -f nvt2.mdp -c nvt1.gro -p topol_cg_sol.top -o nvt2.tpr
$GMX_APP gmx mdrun -ntomp 16 -v -deffnm nvt2

# unrestrained NPT equlibration -- allow the box to relax into the system ~allegro~
$GMX_APP gmx grompp -f npt1.mdp -c nvt2.gro -p topol_cg_sol.top -o npt1.tpr
$GMX_APP gmx mdrun -ntomp 16 -v -deffnm npt1

# NPT production -- assume this is the more interesting ensemble pending some actual lit review
$GMX_APP gmx grompp -f md.mdp -c npt1.gro -p topol_cg_sol.top -o lipase_cg_md.tpr
$GMX_APP gmx mdrun -ntomp 16 -v -deffnm lipase_cg_md -cpi lipase_cg.cpt

