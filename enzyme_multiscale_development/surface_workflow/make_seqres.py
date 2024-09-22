#!/usr/bin/env python3
# copyright 2024 Mason Hooten mason.simulation@gmail.com MIT license

## to prep topol.top suitable for force matching with exclusions
## from topol.top, remove all bonds/angles/dihedrals and add exclusions text

import sys
import re

def get_seqres(pdbfilename):
    '''
    Count atoms in [atoms] directive
    '''
    pdbbasename = pdbfilename.split('.pdb')[0]
    with open(pdbfilename) as f, open(pdbbasename+'_seqres.txt','w') as g:
        for idx,line in enumerate(f,1):
            linetext = line.strip()
#            if linetext and (linetext[0:6] == 'SEQRES'):
            if linetext[0:6] == 'SEQRES':
                g.write(linetext+'\n')
                print(linetext[0:6])
    return
    
if __name__ == '__main__':
    pdbfilename = sys.argv[1]
    get_seqres(pdbfilename)


