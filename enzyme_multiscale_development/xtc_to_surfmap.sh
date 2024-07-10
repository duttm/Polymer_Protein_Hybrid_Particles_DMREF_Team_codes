#!/bin/bash

GMX_CALL="/path/to/containers/gmx21_5.sif gmx"

ATOMSFILE="lipase_CG_atomtypes.dat"
CHARGESFILE="lipase_CG_charges.dat"
SEQRESFILE="lipase_seqres.txt"

MOLNAM="lipase"
MOLSHR="lip"
XTCNAM="lipase_cg_md"
TPRNAM="lipase_cg_md"

for j in {2..10}
do
    MD_DIR=/path/to/data/${MOLNAM}/cg/production/seed_$j

    cp -a $MD_DIR/${XTCNAM}.xtc ./
    cp -a $MD_DIR/${TPRNAM}.tpr ./
    cp -a $MD_DIR/index.ndx ./

    echo 1 1 | $GMX_CALL trjconv -f ${XTCNAM}.xtc -s ${TPRNAM}.tpr -n index.ndx -dt 40000 -b 100000 -e 500000 -fit rot+trans -o ${MOLSHR}cg-seed${j}_100-500ns_dt40ns_fitrottrans.gro -sep

    rm ${XTCNAM}.xtc ${TPRNAM}.tpr index.ndx

    for i in {0..10}
    do
        GROFILE_NOEXT="${MOLSHR}cg-seed${j}_100-500ns_dt40ns_fitrottrans${i}"
        /path/to/run_surfmap.sh \
        $GROFILE_NOEXT $ATOMSFILE $CHARGESFILE $SEQRESFILE
    done
done
