#!/usr/bin/env python3
# make a pdb file out of a gro file
# plus a residue sequence (SEQRES) file .

# SEQRES mimics the section of that name in a standard PDB file.
#TODO SEQRES could be automatically generated from a gro file.

import sys

INFILENAME=sys.argv[1]
INFILENAME_NOEXT=INFILENAME.split('.')[1]
SEQRES_FILE=sys.argv[2]

with open(INFILENAME) as f:
    title=f.readline()
    atom_count=int(f.readline())
    HEADER = f'HEADER    HYDROLASE                               06-DEC-96   XXXX             '
    TITLE = f'TITLE     STRUCTURE OF LIPASE                                                   '
    COMPND = f'COMPND    xx'
    SOURCE = f'SOURCE    xx'
    KEYWDS = f'KEYWDS    xx'
    EXPDTA = f'EXPDTA    X-RAY DIFFRACTION                                                     '
    AUTHOR = f'AUTHOR    xx'
    REVDAT = f'REVDAT    xx'
    JRNL   = f'JRNL      xx'
    
    REMARK1 = f'REMARK   1 xx'
    REMARK2 = f'REMARK   2 xx'
    REMARK3 = f'REMARK   3 xx'
    CRYST1 = f'CRYST1{1.0:>9.3f}{1.0:>9.3f}{1.0:>9.3f}{90.0:>7.2f}{90.0:>7.2f}{90.0:>7.2f} {"P 1":<11}{"1":>4}'
    ORIGX1 = f'ORIGX1      1.000000  0.000000  0.000000        0.00000'
    ORIGX2 = f'ORIGX2      0.000000  1.000000  0.000000        0.00000'
    ORIGX3 = f'ORIGX3      0.000000  0.000000  1.000000        0.00000'
    SCALE1 = f'SCALE1      1.000000  0.000000  0.000000        0.00000'
    SCALE2 = f'SCALE1      0.000000  1.000000  0.000000        0.00000'
    SCALE3 = f'SCALE1      0.000000  0.000000  1.000000        0.00000'
    MODEL = f'MODEL        1'
    ENDMDL = f'ENDMDL'
    MASTER = f'MASTER'
    END = f'END'
    
    print('\n'.join([HEADER, TITLE, COMPND, SOURCE, KEYWDS, EXPDTA, AUTHOR, REVDAT, JRNL, REMARK1, REMARK2, REMARK3]))
    with open(SEQRES_FILE) as g:
        for line in g: print(line.strip())
    print('\n'.join([CRYST1, ORIGX1, ORIGX2, ORIGX3, SCALE1, SCALE2, SCALE3, MODEL]))
    for i in range(atom_count):
        data_line   = f.readline()
        f1 = int(data_line[0:5].strip())        # resnum
        f2 = data_line[5:10].strip()            # restype
        f3 = data_line[10:15].strip()           # atomname
        f4 = int(data_line[15:20].strip())      # atomnum
        f5 = float(data_line[20:28].strip())    # x
        f6 = float(data_line[28:36].strip())    # y
        f7 = float(data_line[36:44].strip())    # z

        record_name="ATOM"
        printline = f'{record_name:<6}{f4:>5} {f3:<4} {f2:<3} A{f1:>4}    {f5*10:>8.3f}{f6*10:>8.3f}{f7*10:>8.3f}{1.00:>6.2f}{0.00:>6.2f}          {f3[0]:>2}'
        print(printline)
    print(f'{"TER":<6}{str(atom_count+1):>5} {"":<4} {f2:<3}  {f1:>4}')

    print('\n'.join([ENDMDL, MASTER, END]))

