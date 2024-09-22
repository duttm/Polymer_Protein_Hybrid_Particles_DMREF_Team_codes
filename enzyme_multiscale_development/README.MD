# Framework for Martini-based Coarse-grained Model of Enzymes

README.TXT
UPDATED 9/21/24 by MH

## Introduction
This repo [a work in progress] contains codes supporting the manuscript noted in the _Citation_ section below. Codes are organized roughly by workflow, including a variety of scripts for file transformation, simulation environment prep, analysis, and plotting.

### Citation 
Hooten M., Murthy N.S., Gormley A., Dutt M. Framework for Martini-based Coarse-grained Model of Enzymes: Development and Experimental Validation. 
[arXiv DOI](https://doi.org/[arxivdoi])

### Simulation
Example scripts and parameter files for all-atom (AA) and coarse-grained (CG) simulation of lipase and dehalogenase.


### Structure Workflow
Calculation of gyration tensor parameters and pair distribution functions.

### Surface Workflow - Evaluation of Simulated Enzyme Surfaces
Production calculation of surface profiles is done using the shell program `xtc_to_surfmap.sh` (or `xtc_to_surfmap-AA.sh`). This program resamples the given trajectories, preprocesses the sample frames, and runs SURFMAP.

Original codes in this repo are:

+ `xtc_tosurfmap.sh`
+ `run_surfmap.sh`
+ `gro2pdbcoords.py`
+ `gro2pqr.py`

#### xtc_to_surfmap.sh
Inputs:

1. `atypes.dat` - File with atom types (one per line) corresponding to the input trajectory. Used by `gro2pqr.py`.
2. `charges.dat` - File with atom charges (one per line) corresponding to the input trajectory. Used by `gro2pqr.py`.
3. `seqres.txt` - File containing only the SEQRES section of the target PDB-formatted structure file.
4. Hard-coded strings related to project naming conventions and directory structure.
5. One hundred (100) simulation trajectories in `.xtc` format containing only the protein atoms of interest, i.e. with solvent and ions excluded.

Outputs:

1. SURFMAP outputs for 100 sampled frames. 100 directories with calculated surface electrostatics, 100 directories with calculated surface values for all other properties.

Dependencies:

1. `run_surfmap.sh`

Pseudocode:

+ for seednum in 1 to 10
    + run `gmx trjconv -f trajectory_seednum.xtc -o input.gro -sep` to sample ten frames
    + for framenum in 1 to 10
        + run `run_surfmap.sh seednum_framenum.gro atypes.dat charges.dat seqres.txt` to produce SURFMAP output

#### run_surfmap.sh
Inputs:

1. `input.gro` - A Gromacs trajectory frame
1. `atypes.dat` 
1. `charges.dat` 
1. `seqres.txt`

Dependencies:

1. `gro2pdbcoords.py`
2. `gro2pqr.py`
3. `SURFMAP` installed in a suitable Python environment. `pip freeze` output for the environment used in this study is available in `requirements_envsurfmap.txt`.

Outputs:

1. SURFMAP outputs for a single trajectory frame

Pseudocode:

1. run `gro2pdbcoords.py` to create a PDB file
1. run `gro2pqr.py` to create a PQR file
3. run `SURFMAP`

#### gro2pdbcoords.py
Inputs:

1. `input.gro`
2. `seqres.txt`

Outputs:

1. `input.pdb`

#### gro2pqr.py
Inputs:

1. `input.gro`
2. `atypes.dat`
3. `charges.dat`

Outputs:

1. `input.pqr`

### Figures
Scripts used to produce figures from the manuscript.

## Acknowledgements
This work was funded by the National Science Foundation (DMREF-2118860, DMR-1654325, OAC-1835449 and CBET-2309852). The authors acknowledge the use of computational resources enabled via an allocation from NSF ACCESS (allocation DMR-140125).

