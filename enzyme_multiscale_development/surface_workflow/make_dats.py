#!/usr/bin/env python3
# copyright 2024 Mason Hooten mason.simulation@gmail.com MIT license

## to prep topol.top suitable for force matching with exclusions
## from topol.top, remove all bonds/angles/dihedrals and add exclusions text

import sys
import re

def get_directivemap(my_filename):
    '''
    Parse Gromacs top file into dict with line numbers for each directive e.g. [atoms]
    Ignores comments (; ...) and excludes special header statements (# ...)
    map = {0:{'name':'string','start':int,'end':int,'linelist':[ints]},1:...}
    '''
    directiveidx = -1
    directivemap = {}
    directivename = ''
    reading_directive = 0
    prev=None
    p = re.compile(r"\[([A-Za-z0-9_\s]+)\]")
    with open(my_filename) as f:
        for idx,line in enumerate(f,1):
            linetext = line.strip()
            if linetext:
                if linetext[0] == '[': # found a directive
                    reading_directive = 1 # in reading mode now
                    directiveidx = directiveidx+1 # increment directive counter
                    if directiveidx > 0 and not directivemap[directiveidx-1]['end']:
                        directivemap[directiveidx-1]['end'] = idx-1 # close previous
                    directivename = p.match(linetext)[1].strip() # get new name
                    directivemap[directiveidx] = {} # init map of the new name
                    directivemap[directiveidx]['name'] = directivename # name
                    directivemap[directiveidx]['start'] = idx # start line
                    directivemap[directiveidx]['end'] = None # start line
                if linetext[0] == ';': continue
                if linetext[0] == '#' and reading_directive:
                    reading_directive = 0
                    if directivemap[directiveidx]['end'] == None:
                        directivemap[directiveidx]['end'] = idx-1 
    if not directivemap[directiveidx]['end']: directivemap[directiveidx]['end'] = idx
    for key in directivemap:
        print(f"{directivemap[key]['name']} {directivemap[key]['start']} {directivemap[key]['end']}")
        directivemap[key]['linelist'] = \
            [ i for i in range(directivemap[key]['start'],directivemap[key]['end']) ]
    return directivemap

def parse_atoms_in_topol(my_filename,my_directivemap):
    '''
    Count atoms in [atoms] directive
    '''
    atomsct = 0
    atoms_lines = []
    atomtypes = []
    atomcharges = []
    for key in my_directivemap:
        if my_directivemap[key]['name'] == 'atoms':
            atoms_lines = atoms_lines + my_directivemap[key]['linelist']
    with open(my_filename) as f:
        for idx,line in enumerate(f,1):
            if idx in atoms_lines:
                linetext = line.strip()
                if linetext and linetext[0] in '0123456789':
                    linedata = linetext.split()
                    atomtypes.append(linedata[1])
                    atomcharges.append(linedata[7])
                    atomsct = atomsct+1
    return atomsct, atomtypes, atomcharges
    
if __name__ == '__main__':
    topol_top = sys.argv[1]
    topol_basename = topol_top.split('.top')[0]
    mymap = get_directivemap(topol_top)
    atomsct,atomtypes,atomcharges = parse_atoms_in_topol(topol_top,mymap)
    with open(topol_basename+'_atomtypes.dat','w') as f:
        for i in atomtypes:
            f.write(i+'\n')
    with open(topol_basename+'_charges.dat','w') as f:
        for i in atomcharges:
            f.write(i+'\n')


