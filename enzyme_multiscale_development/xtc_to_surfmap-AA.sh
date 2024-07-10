#!/bin/bash

GMX_CALL="/path/to/containers/gmx21_5.sif gmx"

ATOMSFILE="lipase_AA_atomtypes.dat"
CHARGESFILE="lipase_AA_charges.dat"
SEQRESFILE="lipase_seqres.txt"

MOLNAM="lipase"
MOLSHR="lip"
XTCNAM="lipaa-protein"
TPRNAM="lipaa-protein"
NDXNAM="lipaa-protein"

for j in {1..10}

do
    MD_DIR=/path/to/data/${MOLNAM}/aa/production/seed_$j

    echo 1 1 | $GMX_CALL trjconv -f ${MD_DIR}/${XTCNAM}.xtc -s ${MD_DIR}/${TPRNAM}.tpr -n ${MD_DIR}/${NDXNAM}.ndx -dt 10000 -b 400000 -e 500000 -fit rot+trans -o ${MOLSHR}aa-seed${j}_400-500ns_dt10ns_fitrottrans.gro -sep

    for i in {0..10}
    do
        GROFILE_NOEXT="${MOLSHR}aa-seed${j}_400-500ns_dt10ns_fitrottrans${i}"
        /path/to/run_surfmap-AA.sh \
        $GROFILE_NOEXT $ATOMSFILE $CHARGESFILE $SEQRESFILE
    done
done
