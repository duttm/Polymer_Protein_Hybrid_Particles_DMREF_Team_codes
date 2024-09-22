#!/usr/bin/env python3

import sys

residues = []
target_residues = ['ARG','ASP','GLU','HIS','LYS']
with open(sys.argv[1]) as file:
  for idx,line in enumerate(file):
    if idx>0:
        data = line.split()
        if data[2] != 'Inf' and data[3] != 'NA':
            residue_strings = data[3:]
            add_residues = [ i.strip(',') for i in residue_strings ]
#            print(add_residues)
            residues = residues + add_residues
        
        residues = list(set(residues))
#        print(residues)
        resnums = sorted([ int(i.split('_')[1]) for i in residues if i.split('_')[0] in target_residues ])
        resnums = [ str(i) for i in resnums ]
        
print(' '.join(resnums))

#for i in sorted(resnums):
#    print(i)
