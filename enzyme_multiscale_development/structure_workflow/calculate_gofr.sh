#!/bin/bash
# copyright 2023 Mason Hooten mason.simulation@gmail.com MIT license

# calculate pair distribution using the A+T code. 

# call this script with $GRO_TRAJ as a command line input.
# GRO_TRAJ should be a gro file with extension i.e. config.gro

# in the code below, pair_distribution.py is a rewrite of 
# Michael P. Allen and Dominic J. Tildesley's pair_distribution.py.


GRO_TRAJ=$1

python3 ~/cod/file_transformation/gro2AT-cnf.py $GRO_TRAJ
mkdir ./gr_$1
mv cnf.* ./gr_$1
cd gr_$1

# This is a call to a modified version of A+T's pair distribution calculation.
# "dr" is the resolution of the distribution.
# "CNF_PRECISION" is a bookkeeping thing. it is the number of places in the
# the filenames that will be created from the frames contained in your .gro 
# trajectory file. A gro file of 10000 frames should used CNF_PRECISION>=5, so 
# that files can be generated like "cnf.00001" thru "cnf.10000"  
echo '{"dr":0.002,"CNF_PRECISION":5}' | python3 ~/cod/analysis/pair_distribution.py

