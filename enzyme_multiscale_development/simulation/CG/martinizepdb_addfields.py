#!/usr/bin/env python3
import sys

INFILENAME=sys.argv[1]
ATOM_COUNT=int(sys.argv[2])
SEQRES_FILE=sys.argv[3]
INFILENAME_NOEXT=INFILENAME.split('.')[1]

with open(INFILENAME) as f:
    title=f.readline()
    HEADER = f'HEADER {INFILENAME_NOEXT}.pdb from gro using gro2pdb'
    TITLE = f'TITLE'
    COMPND = f'COMPND'
    SOURCE = f'SOURCE'
    KEYWDS = f'KEYWDS'
    EXPDTA = f'EXPDTA'
    AUTHOR = f'AUTHOR'
    REVDAT = f'REVDAT   1   DD-MMM-YY xPDB    0'
    JRNL   = f'JRNL'
    REMARK1 = f'REMARK   1'
    REMARK2 = f'REMARK   2'
    REMARK3 = f'REMARK   3'
    CRYST1 = f'CRYST1'
    ORIGX1 = f'ORIGX1      1.000000  0.000000  0.000000        0.00000'
    ORIGX2 = f'ORIGX2      0.000000  1.000000  0.000000        0.00000'
    ORIGX3 = f'ORIGX3      0.000000  0.000000  1.000000        0.00000'
    SCALE1 = f'SCALE1'
    SCALE2 = f'SCALE2'
    SCALE3 = f'SCALE3'
    MODEL  = f'MODEL   1'
    ENDMDL = f'ENDMDL'
    MASTER = f'MASTER'
    END = f'END'
    print('\n'.join([HEADER, TITLE, COMPND, SOURCE, KEYWDS, EXPDTA, AUTHOR, REVDAT, REMARK1, REMARK2, REMARK3]))
    with open(SEQRES_FILE) as g:
        print(''.join(g.readlines()))
    print('\n'.join([CRYST1, ORIGX1, ORIGX2, ORIGX3, SCALE1, SCALE2, SCALE3, MODEL]))
    for i in range(ATOM_COUNT):
        data_line   = f.readline().strip()

#        f1 = int(data_line[0:5].strip())        # resnum
#        f2 = data_line[5:10].strip()            # restype
#        f3 = data_line[10:15].strip()           # atomname
#        f4 = int(data_line[15:20].strip())      # atomnum
#        f5 = float(data_line[20:28].strip())    # x
#        f6 = float(data_line[28:36].strip())    # y
#        f7 = float(data_line[36:44].strip())    # z
#        one_atom = [f1,f2,f3,f4,f5,f6,f7]
#        
#        record_name="ATOM"
#        printline = f'{record_name:<6}{f4:>5} {f3:<4} {f2:<3}  {f1:>4}    {f5*10:>8.3f}{f6*10:>8.3f}{f7*10:>8.3f}{1.00:>6.2f}{0.00:>6.2f}'
#        print(printline)
        print(data_line)
    print('\n'.join([ENDMDL,MASTER, END]))

